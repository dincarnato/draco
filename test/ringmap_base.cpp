#include <ringmap_data.hpp>
#include <ringmap_matrix.hpp>
#include <tokenizer_iterator.hpp>

#include <fstream>
#include <string>
#include <cassert>

#include <armadillo>

static constexpr const char* sequence = "ACCACCAACA";

arma::Mat<unsigned>
readRingmapFile(const char* filename)
{
    unsigned begin = std::numeric_limits<unsigned>::max();
    unsigned end = 0;
    std::ifstream file(filename);
    unsigned nReads = 0;
    for(std::string line; std::getline(file, line); ++nReads)
    {
        auto&& tokenizer = makeTokenizerIterator(line, '\t');
        begin = std::min(begin, tokenizer++.tou() - 1);
        end = std::max(end, tokenizer.tou());
    }

    file.clear();
    file.seekg(0);

    arma::Mat<unsigned> mat(nReads, end - begin, arma::fill::zeros);
    unsigned readIndex = 0;
    for(std::string line; std::getline(file, line); ++readIndex)
    {
        auto&& tokenizer = makeTokenizerIterator(line, '\t');
        unsigned currentBegin = tokenizer.tou() - 1;
        tokenizer += 2;
        unsigned baseIndex = 0;
        for(char value : *tokenizer)
        {
            if(value == '1')
                mat(readIndex, baseIndex + currentBegin - begin) = 1;
            ++baseIndex;
        }
    }

    return mat;
}

int main(int argc, const char* argv[])
{
    assert(argc == 2);
    const char* const filename = argv[1];

    Args args;
    RingmapData ringmapData(filename, sequence, args);
    const RingmapMatrix& matrix = ringmapData.data();
    auto armaMatrix = readRingmapFile(filename);

    assert(matrix.rows_size() == armaMatrix.n_rows);
    assert(matrix.cols_size() == armaMatrix.n_cols);
    assert(matrix == matrix);
    assert(not(matrix != matrix));

    /* Element-wise check */
    for(unsigned row = 0; row < matrix.rows_size(); ++row)
    {
        for(unsigned col = 0; col < matrix.cols_size(); ++col)
            assert(matrix(row, col) == armaMatrix(row, col));
    }

    /* Row-wise check */
    {
        unsigned rowIndex = 0;
        for(auto&& row : matrix.rows())
        {
            assert(row.sum() == arma::sum(armaMatrix.row(rowIndex)));
            assert(row.mean() == arma::mean(arma::conv_to<arma::mat>::from(armaMatrix).row(rowIndex)));

            unsigned colIndex = 0;
            for(bool mutated : row)
                assert(mutated == armaMatrix(rowIndex, colIndex++));
            ++rowIndex;
        }
    }

    /* Col-wise check */
    {
        unsigned colIndex = 0;
        for(auto&& col : matrix.cols())
        {
            assert(col.sum() == arma::sum(armaMatrix.col(colIndex)));
            assert(col.mean() == arma::mean(arma::conv_to<arma::mat>::from(armaMatrix).col(colIndex)));

            unsigned rowIndex = 0;
            for(bool mutated : col)
                assert(mutated == armaMatrix(rowIndex++, colIndex));
            ++colIndex;
        }
    }

    /* Transposed matrix element-wise check */
    {
        decltype(armaMatrix) armaTransposed = armaMatrix.t();
        auto transposed = matrix.t();

        assert(transposed.rows_size() == armaTransposed.n_rows);
        assert(transposed.cols_size() == armaTransposed.n_cols);

        for(unsigned row = 0; row < transposed.rows_size(); ++row)
        {
            for(unsigned col = 0; col < transposed.cols_size(); ++col)
                assert(transposed(row, col) == armaTransposed(row, col));
        }
    }
}
