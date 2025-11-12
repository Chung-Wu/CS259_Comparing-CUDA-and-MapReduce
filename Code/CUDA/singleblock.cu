#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

// ===================== Macro for CUDA error checking =====================
#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s %s %d\n",
                cudaGetErrorString(code), file, line);
        exit(code);
    }
}

// ===================== GPU Kernel =====================
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
                default:  val = 4; break;  
            }
            code = code * 5 + val;  // base-5 
        }
        atomicAdd(&counts[code], 1);
    }
}

// ===================== Host main function =====================
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

    //  initialize GPU Memory
    int num_kmers = 1;
    for (int i = 0; i < k; i++) num_kmers *= 5;  // 5^k

    int *d_counts;
    CUDA_CHECK(cudaMalloc(&d_counts, num_kmers * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_counts, 0, num_kmers * sizeof(int)));

    const size_t CHUNK_SIZE = 50 * 1024 * 1024; // 50 MB each time
    char *d_seq;
    CUDA_CHECK(cudaMalloc(&d_seq, CHUNK_SIZE));

    int threadsPerBlock = 256;
    size_t processed = 0;
    
    // start measuring the time
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    // 4. processing different chunk of data
    while (processed < (size_t)n) {
        size_t len = min(CHUNK_SIZE, (size_t)n - processed);
        CUDA_CHECK(cudaMemcpy(d_seq, data.c_str() + processed, len, cudaMemcpyHostToDevice));

        int numBlocks = 1;

        kmer_count_kernel<<<numBlocks, threadsPerBlock>>>(d_seq, len, k, d_counts);
        CUDA_CHECK(cudaDeviceSynchronize());

        processed += CHUNK_SIZE - (k - 1);
    }

    // 5. stop timing
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "GPU processing time: " << milliseconds / 1000.0 << " sec" << endl;

    int* counts = new int[num_kmers];
    CUDA_CHECK(cudaMemcpy(counts, d_counts, num_kmers * sizeof(int), cudaMemcpyDeviceToHost));

    cout << "\n=== Non-zero k-mer counts ===" << endl;
    for (int i = 0; i < num_kmers; i++) {
        if (counts[i] > 0)
            printf("k-mer code %3d : %d\n", i, counts[i]);
    }

    delete[] counts;
    cudaFree(d_seq);
    cudaFree(d_counts);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cout << "\nDone.\n";
    return 0;
}
