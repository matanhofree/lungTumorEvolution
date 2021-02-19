/*
 * fastTableRead.cpp
 *
 *  Created on: Aug 5th, 2014
 *      Author: mhofree
 *
 */
#include "mex.h"
#include <iterator>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sys/stat.h>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
//#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/detail/ios.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/filesystem.hpp>
#include <cassert>

// TODO: fix empty delim bug!
#define xMM(A,m,n,D) (A[m+(D*n)])

using namespace std;

//template<typename T>
//class Allocator {
//public :
//    //    typedefs
//    typedef T value_type;
//    typedef value_type* pointer;
//    typedef const value_type* const_pointer;
//    typedef value_type& reference;
//    typedef const value_type& const_reference;
//    typedef std::size_t size_type;
//    typedef std::ptrdiff_t difference_type;
//
//public :
//    //    convert an allocator<T> to allocator<U>
//    template<typename U>
//    struct rebind {
//        typedef Allocator<U> other;
//    };
//
//public :
//    inline explicit Allocator() {}
//    inline ~Allocator() {}
//    inline explicit Allocator(Allocator const&) {}
//    template<typename U>
//    inline explicit Allocator(Allocator<U> const&) {}
//
//    //    address
//    inline pointer address(reference r) { return &r; }
//    inline const_pointer address(const_reference r) { return &r; }
//
//    //    memory allocation
//    inline pointer allocate(size_type cnt,
//                            typename std::allocator<void>::const_pointer = 0) {
//        //return reinterpret_cast<pointer>(::operator mxMalloc(cnt * sizeof (T)))
//        return reinterpret_cast<pointer>(mxMalloc(cnt * sizeof (T)));
//    }
//    inline void deallocate(pointer p, size_type) {
//        mxFree(p);
//    }
//
//    //    size
//    inline size_type max_size() const {
//        return std::numeric_limits<size_type>::max() / sizeof(T);
//    }
//
//    //    construction/destruction
//    inline void construct(pointer p, const T& t)
//    {
//        mxMalloc(p);
//        *p = T(t);
//    }
//    inline void destroy(pointer p) { p->~T(); }
//
//    inline bool operator==(Allocator const&) { return true; }
//    inline bool operator!=(Allocator const& a) { return !operator==(a); }
//};    //    end of class Allocator
//struct matlab_alloc : Allocator<char> { };

uint32_T countLines(const char * fname)
{

	ifstream myfile(fname);
	size_t numLines = 0;

	if (myfile.is_open() && myfile.good())
	{
		// new lines will be skipped unless we stop it from happening:
		myfile.unsetf(std::ios_base::skipws);

		// count the newlines with an algorithm specialized for counting:
		numLines = std::count(
				std::istream_iterator<char>(myfile),
				std::istream_iterator<char>(),'\n');
	}
    return numLines;
}


int readFileData(char* fileName,int &rowLen, char* delimChars, int skipLines, char* commentChar, vector<string> &errorLines,vector<vector<string> > &dataTable)
{
    // Read gzipped files not yet supported
    uint32_T cntLines = 0;
    boost::iostreams::filtering_istream inStreamFilter;

    std::string::size_type idx;
    idx = string(fileName).rfind('.');

    bool isGzip = false;
    if(idx != std::string::npos)
    {
        std::string extension = string(fileName).substr(idx+1);
        // mexPrintf("Ext:(%s)\n",extension.c_str());
        // isGzip = extension.compare("gz") == 0;
        isGzip = strncmp(extension.c_str(),"gz",2) == 0;

        if (isGzip)
        {
            mexPrintf("Gzip not supported\n");
            return -1;
        }
    }


    std::ios_base::openmode fopenMode = std::ios_base::in;
    if (isGzip)
    {
        fopenMode = fopenMode | std::ios_base::binary;
    }

    boost::iostreams::file_source inF(fileName, fopenMode);
    if (!inF.is_open())
    {
        return -2;
    }
//    if (isGzip) {
//        mexPrintf("Pre-gzip decomp\n");
//        inStreamFilter.push(boost::iostreams::gzip_decompressor<(std::allocator)Allocator<char> >());
//        mexPrintf("gzip decomp\n");
//    }
    inStreamFilter.push(inF);


    //istream inStream(&inF);
    //std::istream inStream(&inStreamFilter);
    string line;
    while (inStreamFilter) {

        getline(inStreamFilter,line);
        if ( inStreamFilter.eof() ) break;

        //mexPrintf("Txt: %s\n",line.c_str());
        //mexPrintf("Line: %d\n",cntLines);

        vector<string> tokens;
        boost::split(tokens, line, boost::is_any_of(delimChars));
// Utils::Tokenize(line, tokens, " \t");
// Tokenize
//					boost::tokenizer<boost::char_separator<char>> tok(line, sep);
//					for (boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg){
//						tokens.push_back()
//					}

        // Skip header lines
        if (skipLines > 0)
        {
            errorLines.push_back(line);
            --skipLines;
            continue;
        }

        bool foundComment = 0;
        if (commentChar)
        {
            // Simple first in line comment masker
            if (line[0] == commentChar[0])
            {
                errorLines.push_back(line);
                continue;
            }


//            mexPrintf("Comment char:%s\n",commentChar);
//            for (vector<string>::iterator dataItCol = tokens.begin(); dataItCol != tokens.end(); ++dataItCol)
//            {
//                if (dataItCol->find(commentChar) != std::string::npos)
//                {
//                    errorLines.push_back(tokens.
//                    tokens.erase(dataItCol,tokens.end());
//                    foundComment = 1;
//                    errorLines.push_back(line);
//                    break;
//                }
//            }
//            if ( foundComment && rowLen == 0 )
//            {
//                continue;
//            }
        }

        if ( rowLen == 0 )
        {
            rowLen = tokens.size();
//            mexPrintf("Line #: %d, skipped %d\n",cntLines,skipLines);
//            mexPrintf("Txt:%s\n",line.c_str());

        };

        if (rowLen == tokens.size())
        {
            dataTable.push_back(tokens);
            ++cntLines;
        }
        else
        {
            errorLines.push_back(line);
        }
    }
    return cntLines;
}

//int ends_with(const char* name, const char* extension, size_t length)
//{
//    const char* ldot = strrchr(name, '.');
//    if (ldot != NULL)
//    {
//        if (length == 0)
//            length = strlen(extension);
//        return strncmp(ldot + 1, extension, length) == 0;
//    }
//    return 0;
//}


/*******************************************************************************
 **entry point for the mex call
 **nlhs - number of outputs
 **plhs - pointer to array of outputs
 **nrhs - number of inputs
 **prhs - pointer to array of inputs
 *******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*****************************************************************************
	 * A very simple reader function for reading a text table from file
	 * Input:
	 * flag: (1) Read text output will be a simple cell array
	 *       (2) (TODO) Read data frame
	 *       (3)
	 * filename: name of file to write
	 * colHeaderLines:
	 * rowHeadersLines:
	 * delim:
	 *
	 * Output:
	 * outCellData
	 * outErrorLines
	 *
	 * TODO:
	 * - add simple error checking.
	 *****************************************************************************/

	const mxArray* mxFlag;
	const mxArray* mxFileName;
	const mxArray* mxData = NULL;
	const mxArray* mxDataRow;
	const mxArray* mxDataCol;
	const mxArray* mxTxtRow;
	const mxArray* mxTxtCol;
	const mxArray* mxRowNum = NULL;
	const mxArray* mxDelim = NULL;
    const mxArray* mxCommentChar = NULL;


    int skipLines = 0;
    int rowNames = 0;
    int colNames = 0;
    char* commentChar = "";
	switch (nrhs)
	{
		case 2:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			/*printf("Must be H(X), calculateProbability(X), merge(X), normaliseArray(X)\n");*/
			break;
		}
		case 3:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			mxRowNum = prhs[2];
			/*printf("Must be H(X), calculateProbability(X), merge(X), normaliseArray(X)\n");*/
			break;
		}
		case 4:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			mxRowNum = prhs[2];
			mxDelim = prhs[3];
			break;
		}
		case 5:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			mxRowNum = prhs[2];
			mxDelim = prhs[3];
            if (prhs[4]) skipLines = *mxGetPr(prhs[4]);
			break;
		}
		case 7:
		{
            mxFlag = prhs[0];
            mxFileName = prhs[1];
            mxRowNum = prhs[2];
            mxDelim = prhs[3];
            if (prhs[4]) skipLines = *mxGetPr(prhs[4]);
            if (prhs[5]) rowNames = *mxGetPr(prhs[5]);
            if (prhs[6]) colNames = *mxGetPr(prhs[6]);
            break;
		}
        case 8:
        {
            mxFlag = prhs[0];
            mxFileName = prhs[1];
            mxRowNum = prhs[2];
            mxDelim = prhs[3];
            if (prhs[4]) skipLines = *mxGetPr(prhs[4]);
            if (prhs[5]) rowNames = *mxGetPr(prhs[5]);
            if (prhs[6]) colNames = *mxGetPr(prhs[6]);
            if (prhs[7]) commentChar = mxArrayToString(prhs[7]);
            break;
        }
		default:
		{
			mexPrintf("Error incorrect number of arguments, format is SimpleByteReadWrite(fcn,data)\n");
			return;
		}
	}

	int flag = *mxGetPr(mxFlag);

	switch (flag)
	{
		case 1: /* read simple text data */
		{

			char* fileName = mxArrayToString(mxFileName);
            int rowLen = 0;
            vector<string> errorLines;
            vector<vector<string> > dataTable;

            char * delimChars;
            if (mxDelim != NULL)
            {
                delimChars = mxArrayToString(mxDelim);
                // mexPrintf("Delim:(%s)\n",delimChars);
            }

            int errorState = readFileData(fileName,rowLen,delimChars,skipLines,commentChar,errorLines,dataTable);
            if ( errorState < 0 )
            {
                mexPrintf("Error no data read\n");
            }
            mxArray* outCellData = mxCreateCellMatrix(dataTable.size(),rowLen);
            mxArray* outErrorLines = mxCreateCellMatrix(errorLines.size(),1);
            mwIndex cellPos[2] = {0,0};

            mwIndex cntCopy = 0;
            for (vector<vector<string> >::iterator dataItRow = dataTable.begin(); dataItRow != dataTable.end();++dataItRow)
            {
                cellPos[1]=0;
                for (vector<string>::iterator dataItCol = dataItRow->begin(); dataItCol != dataItRow->end();++dataItCol)
                {
                    const char * zstr = dataItCol->c_str();
                    cntCopy = mxCalcSingleSubscript(outCellData, rowLen, cellPos);

                    mxSetCell(outCellData, cntCopy, mxCreateString(zstr));
                    ++cellPos[1];
                }
                ++cellPos[0];
            }

            cntCopy = 0;
            for (vector<string>::iterator dataItCol = errorLines.begin(); dataItCol != errorLines.end();++dataItCol)
            {
                const char * zstr = dataItCol->c_str();
                mxSetCell(outErrorLines, cntCopy++, mxCreateString(zstr));
            }


            plhs[0] = outCellData;
            plhs[1] = outErrorLines;
            break;
		}
		case 2: /* read simple numeric only data with top and bottom headers */
		{
            char* fileName = mxArrayToString(mxFileName);
            int dataDim = 0;
            vector<string> errorLines;
            vector<vector<string> > dataTable;

            char * delimChars;
            if (mxDelim != NULL)
            {
                delimChars = mxArrayToString(mxDelim);
                // mexPrintf("Delim:(%s)\n",delimChars);
            }

            int errorState = readFileData(fileName,dataDim,delimChars,skipLines,commentChar,errorLines,dataTable);

            mwIndex nDim = dataDim - rowNames;
            mwIndex nRow = dataTable.size() - colNames;

            mxArray* outData = mxCreateDoubleMatrix(nRow,nDim,mxREAL);;
            mxArray* outErrorLines = mxCreateCellMatrix(errorLines.size(),1);
            mxArray* outColNames;
            mxArray* outRowNames;
            if (colNames>0)
            {
                outColNames = mxCreateCellMatrix(nDim,colNames);
            }
            if (rowNames>0)
            {
                outRowNames = mxCreateCellMatrix(nRow,rowNames);
            }


            double nan = mxGetNaN();
            double* outDataPr = mxGetPr(outData);

            mwIndex cellPos[2] = {0,0};
            mwIndex dataPos[2] = {0,0};
            mwIndex nColHeader = 0;
            mwIndex nRowHeader = 0;
            mwIndex cntCopy = 0;



            for (vector<vector<string> >::iterator dataItRow = dataTable.begin(); dataItRow != dataTable.end();++dataItRow)
            {
                cellPos[1]=0;
                dataPos[1]=0;

                for (vector<string>::iterator dataItCol = dataItRow->begin(); dataItCol != dataItRow->end();++dataItCol)
                {

                    if (rowNames > 0 && cellPos[1] < rowNames)
                    {
                        // Row header
                        if (cellPos[0] >= colNames ) {
                            const char *zstr = dataItCol->c_str();
                            mxSetCell(outRowNames, nRowHeader, mxCreateString(zstr));
                            assert(nRowHeader < nRow);
                            ++nRowHeader;
                        }
                    }
                    else if (colNames > 0 && cellPos[0] < colNames)
                    {
                        // Col header
                        if (cellPos[1] >= rowNames) {
                            const char *zstr = dataItCol->c_str();
                            mxSetCell(outColNames, nColHeader, mxCreateString(zstr));
                            assert(nColHeader < nDim);
                            ++nColHeader;
                        }
                    }
                    else {
                        cntCopy = mxCalcSingleSubscript(outData, 2, dataPos);

                        try {
                            *(outDataPr + cntCopy) = boost::lexical_cast<double>(*dataItCol);
                            ++dataPos[1];
                        }
                        catch (const boost::bad_lexical_cast &)
                        {
                            *(outDataPr + cntCopy) = nan;
                            ++dataPos[1];
                        }
                    }
                    ++cellPos[1];
                }

                if (cellPos[0] >= colNames) ++dataPos[0];
                ++cellPos[0];
            }

            cntCopy = 0;
            for (vector<string>::iterator dataItCol = errorLines.begin(); dataItCol != errorLines.end();++dataItCol)
            {
                const char * zstr = dataItCol->c_str();
                mxSetCell(outErrorLines, cntCopy++, mxCreateString(zstr));
            }


            plhs[0] = outData;
            plhs[1] = outErrorLines;
            if (rowNames>0) plhs[2] = outRowNames;
            if (colNames>0) plhs[3] = outColNames;

            break;
		}
		default:
		{
			printf("Unrecognised flag %d\n",flag);
			break;
		}/*default*/
	}

	return;
}/*mexFunction()*/




