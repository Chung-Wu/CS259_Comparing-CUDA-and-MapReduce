
#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

// ===================== macro for error checking =====================
#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s %s %d\n",
                cudaGetErrorString(code), file, line);
        exit(code);
    }
}

// ===================== Kernel：k-mer counting =====================
__global__ void kmer_count_kernel(const char* seq, int n, int k, int* counts) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i + k <= n) {
        int code = 0;
        for (int j = 0; j < k; j++) {
            char c = seq[i + j];
            int val;
            switch (c) {
                case 'A': val = 0; break;
                case 'C': val = 1; break;
                case 'G': val = 2; break;
                case 'T': val = 3; break;
                case 'N': val = 4; break;
                default:  val = 4; break; // 其他符號當 N
            }
            code = code * 5 + val; // base-5 編碼
        }
        atomicAdd(&counts[code], 1);
    }
}

int main() {
    ifstream file("../../Data/hg19.fa");
    if (!file.is_open()) {
        cerr << "Error: cannot open file.\n";
        return 1;
    }

    string line, data;
    while (getline(file, line)) {
        if (!line.empty() && line[0] != '>') {
            data += line;
        }
    }
    file.close();

    unsigned long n = data.size();
    int k = 3;
    cout << "Total sequence length: " << n << " bp" << endl;
    
    // allocate GPU ram
    int num_kmers = 1;
    for (int i = 0; i < k; i++) num_kmers *= 5;
    int *d_counts;
    CUDA_CHECK(cudaMalloc(&d_counts, num_kmers * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_counts, 0, num_kmers * sizeof(int)));

    const size_t CHUNK_SIZE = 50 * 1024 * 1024; // 50MB
    const int NUM_STREAMS = 2;

    char* h_seq[NUM_STREAMS];
    char* d_seq[NUM_STREAMS];
    cudaStream_t streams[NUM_STREAMS];

    for (int i = 0; i < NUM_STREAMS; i++) {
        CUDA_CHECK(cudaHostAlloc(&h_seq[i], CHUNK_SIZE, cudaHostAllocDefault)); // pinned host mem
        CUDA_CHECK(cudaMalloc(&d_seq[i], CHUNK_SIZE)); // GPU buffer
        CUDA_CHECK(cudaStreamCreate(&streams[i]));
    }

    int threadsPerBlock = 256;
    size_t processed = 0;
    int stream_id = 0;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    while (processed < (size_t)n) {
        size_t len = min(CHUNK_SIZE, (size_t)n - processed);

        memcpy(h_seq[stream_id], data.c_str() + processed, len);

        CUDA_CHECK(cudaMemcpyAsync(d_seq[stream_id], h_seq[stream_id],
                                   len, cudaMemcpyHostToDevice, streams[stream_id]));

        int numBlocks = 1;

        kmer_count_kernel<<<numBlocks, threadsPerBlock, 0, streams[stream_id]>>>(
            d_seq[stream_id], len, k, d_counts);

        processed += CHUNK_SIZE - (k - 1);
        stream_id = (stream_id + 1) % NUM_STREAMS;
    }

    CUDA_CHECK(cudaDeviceSynchronize());

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "GPU pipelined (pinned+streams) processing time: "
         << milliseconds / 1000.0 << " sec" << endl;

    int* counts = new int[num_kmers];
    CUDA_CHECK(cudaMemcpy(counts, d_counts, num_kmers * sizeof(int), cudaMemcpyDeviceToHost));

    cout << "\n=== Non-zero k-mer counts ===" << endl;
    for (int i = 0; i < num_kmers; i++) {
        if (counts[i] > 0)
            printf("k-mer code %3d : %d\n", i, counts[i]);
    }

    delete[] counts;
    for (int i = 0; i < NUM_STREAMS; i++) {
        cudaFreeHost(h_seq[i]);
        cudaFree(d_seq[i]);
        cudaStreamDestroy(streams[i]);
    }
    cudaFree(d_counts);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cout << "\nDone.\n";
    return 0;
}
