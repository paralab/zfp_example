#include <iostream>
#include <cmath>
#include <iomanip>
#include "zfp.h"

using namespace std;

int main() {

    // ***************************************
    // Initialization
    // ***************************************

    // define and initialize an array
    int nx = 4, ny = 7;
    double a[ny][nx];

    cout << "\noriginal array:" << endl;
    for (int j = 0; j < ny; j++){
        for (int i = 0; i < nx; i++) {
            double x = 2.0 * i / nx;
            double y = 3.0 * j / ny;
            a[j][i] = exp(-(x * x + y * y));
            cout << setprecision(16) << a[j][i] << endl;
        }
    }

    printf("\n===============================\n");
    printf("\noriginal array size = %lu bytes\n", sizeof(a));

    // ====================================================================================
    // Compression modes (from the documentation)
    // ====================================================================================

    // 1- Expert Mode
    // 2- Fixed-Rate Mode
    // 3- Fixed-Precision Mode
    // 4- Fixed-Accuracy Mode
    // 5- Reversible Mode

    // 1- Expert Mode: there are 4 parameters that can be set directly: minbits, maxbits, maxprec, minexp. Check
    // Section 6.1 of the documentation for details.

    // 2- Fixed-Rate Mode: It is needed to support random access to blocks, and also is the mode used in the
    // implementation of zfp’s compressed arrays. Fixed-rate mode also ensures a predictable memory/storage footprint,
    // but usually results in far worse accuracy per bit than the variable-rate fixed-precision and fixed-accuracy modes.

    // 3- Fixed-Precision Mode: It is preferable when relative rather than absolute errors matter.

    // 4- Fixed-Accuracy Mode: It gives the highest quality (in terms of absolute error) for a given compression rate,
    // and is preferable when random access is not needed.

    // 5- Reversible Mode: It is the lossless compression method. As with the other compression modes, each block is
    // compressed and decompressed independently, but reversible mode uses a different compression algorithm that
    // ensures a bit-for-bit identical reconstruction of integer and floating-point data.

    // ====================================================================================
    // Syntax for different modes:
    // ====================================================================================

    // 1- Expert Mode
    // int zfp_stream_set_params(zfp_stream* stream, uint minbits, uint maxbits, uint maxprec, int minexp)

    // 2- Fixed-Rate Mode
    // double zfp_stream_set_rate(zfp_stream* stream, double rate, zfp_type type, uint dims, int wra)
    // "wra" should be nonzero for random access.

    // 3- Fixed-Precision Mode
    // uint zfp_stream_set_precision(zfp_stream* stream, uint precision)

    // 4- Fixed-Accuracy Mode
    // double zfp_stream_set_accuracy(zfp_stream* stream, double tolerance)

    // 5- Reversible Mode
    // void zfp_stream_set_reversible(zfp_stream* stream)

    // ====================================================================================
    // Compression
    // ====================================================================================

    // The zfp_field parameter object holds information about the uncompressed array.
    // allocate metadata for the 2D array a[ny][nx]
    uint dims = 2;
    zfp_type type = zfp_type_double;
    zfp_field* field = zfp_field_2d(&a[0][0], type, nx, ny);

    // To specify the compressed array, a zfp_stream object must be allocated.
    // allocate metadata for a compressed stream
    zfp_stream* zfp = zfp_stream_open(nullptr);

    // set compression mode and parameters
    double tolerance = 1e-6;
    double rate = 16;

    // uncomment one of the following compression modes
//    zfp_stream_set_rate(zfp, rate, type, dims, 0);
//    zfp_stream_set_precision(zfp, rate);
    zfp_stream_set_accuracy(zfp, tolerance);
//    zfp_stream_set_reversible(zfp);

    // The compression parameters have now been specified, but before compression can occur a buffer large enough to
    // hold the compressed bit stream must be allocated. Another utility function exists for estimating how many bytes
    // are needed.
    // allocate buffer for compressed data
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    auto buffer = new uchar[bufsize];

//    printf("\nmaximum compression buffer size = %lu bytes\n", bufsize);

    // associate bit stream with allocated buffer
    bitstream* stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);

    // compress entire array
    size_t csize = zfp_compress(zfp, field);

    if (!csize) {
        fprintf(stderr, "compression failed\n");
    }else{
        printf("\ncompression completed!\n");
        printf("compressed array size = %lu bytes\n", csize);
    }

    // ====================================================================================
    // Decompression
    // ====================================================================================

    // rewind compressed stream and decompress array
    zfp_stream_rewind(zfp);
    size_t desize = zfp_decompress(zfp, field);

    if (!desize) {
        fprintf(stderr, "decompression failed\n");
    }else{
        printf("\ndecompression completed!\n");
    }

    printf("\n===============================\n");
    cout << "\ndecompressed array: " << endl;
    for (int j = 0; j < ny; j++){
        for (int i = 0; i < nx; i++) {
            cout << setprecision(16) << a[j][i] << endl;
        }
    }

    // ====================================================================================
    // Clean up and finalize
    // ====================================================================================

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    delete []buffer;

    return 0;
}
