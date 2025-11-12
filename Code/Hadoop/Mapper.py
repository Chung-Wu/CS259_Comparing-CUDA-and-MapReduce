import sys

# def mapper(size: int = 3):
#     for data in sys.stdin:
#         if data.startswith(">"):
#             continue
#         sliding_window(data, size)
#         # brutal(data, size)
#
#
# def sliding_window(data: str, size: int):
#     # data = data.strip()                             # ensure removing the space on both sides
#     print("sidling")
#     data = data.replace("\n", "").strip()
#     if len(data) < size:
#         return
#     window = data[:size]  # first to assign the data from index 0 to index kmer_size
#     print(f"{window}\t1")  # ensure {key, value} pair (contains no extra space)
#     for i in range(1, len(data) - size + 1):
#         # avoid double loop. Instead, since we're merely poping front digit and pushing another digit in the tail,
#         # we can keep the majority of list.
#         # This method can significantly reduce the time complexity.
#
#         window = window[1:] + data[i + size - 1]
#         print(f"{window}\t1")
#
#
# def brutal(data: str, size: int):
#     print("brutal")
#     for i in range(len(data) - size):
#         # from i to (lenth of data - kmer_size)
#         for j in range(
#             size
#         ):  # here is for second loop, changing from 0 to lengt of kmer_size
#             print(data[i + j], end="")
#         print("\t1")  # ensure the {key, value} pair
#
#
# if __name__ == "__main__":
#     try:
#         mapper(int(sys.argv[1]))  # To set the size for kmer
#     except:
#         mapper()  # the size is default 3


def mapper(size: int = 3):
    buffer = ""                                 # for buffering, which is used to store the previous input tail.
    for data in sys.stdin:
        if data.startswith(">"):                # skip the first line, which begins with >gi...
            continue
        data = data.strip()                     # for removing the head and tail spaces
        if not data:                            # if there is nothing after removing these spaces, we pass this input (just for double check)
            continue

        data = buffer + data                    # concatenate buffer (previous ending before '\n') and data (current input beginning)

        if len(data) >= size:                   # if the length is larger than kmer_size, then do the text seperating
            sliding_window(data, size)

        buffer = data[-(size - 1) :]            # update buffer for next input after '\n'

def sliding_window(data: str, size: int):
    for i in range(len(data) - size + 1):
        print(f"{data[i : i + size]}\t1")


if __name__ == "__main__":
    try:
        mapper(int(sys.argv[1]))  # To set the size for kmer
    except:
        mapper()  # default k-mer size = 3
