#include "matrix_io.h"

#include <regex>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cmath>
#include <iomanip>
#include <iterator>

#include "../../helpers.hpp"

void print_matrix_label(std::ostream& ostr, const Matrix& matrix)
{
	switch(matrix.type) {
		case DYNAMIC: 
			{
				ostr << "dynamic";

				break;
			}
		case F:
			{
				auto m = dynamic_cast<const FMatrix&>(matrix);	
				ostr << "F_{" << m.k << "," << m.b << "}";
				break;
			}
	    case DIAGONAL:
            {
                ostr << "diagonal";

                break;
            }
		case ONE_B:
			{
				auto m = dynamic_cast<const OneBMatrix&>(matrix);	
				ostr << "1_" << m.b;
				break;
			}
		case TERMINATOR:
			{
				ostr << "T";
				break;
			}
		case ZERO:
			{
				ostr << "0";
				break;
			}
	}
}

void write_matrix(std::ostream& ostr, const Matrix& matrix)
{
	switch(matrix.type) {
		case DYNAMIC: 
			{
				auto m = dynamic_cast<const DynamicMatrix&>(matrix);	
				ostr << "dynamic " << m.n << " " << m.m << " [";
				unsigned long long i_max_row = 1;
				unsigned long long i_max_col = 1;
				i_max_col = i_max_col << matrix.n;
				i_max_row = i_max_row << matrix.m;

				for(unsigned long long i_row = 0; i_row < i_max_row; i_row++)
				{
					ostr << ((i_row == 0) ? "" : ";");
					for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
						ostr << ((i_col == 0) ? "" : ",") << m.get(i_row, i_col);
				}
				ostr << "]";

				break;
			}
		case F:
			{
				auto m = dynamic_cast<const FMatrix&>(matrix);	
				ostr << "F_{" << m.k << "," << m.b << "}";
				break;
			}
        case DIAGONAL:
            {
                auto m = dynamic_cast<const DiagonalMatrix&>(matrix);
                ostr << "diagonal " << m.k << " [";
                unsigned long long i_max_col = 1;
                i_max_col = i_max_col << matrix.n;

                for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
                    ostr << ((i_col == 0) ? "" : ",") << m.get(i_col, i_col);
                ostr << "]";

                break;
            }
		case ONE_B:
			{
				auto m = dynamic_cast<const OneBMatrix&>(matrix);	
				ostr << "1_" << m.b;
				break;
			}
		case TERMINATOR:
			{
				ostr << "T";
				break;
			}
		case ZERO:
			{
				auto m = dynamic_cast<const ZeroMatrix&>(matrix);	
				ostr << "0_{" << m.n << "," << m.m << "}";
				break;
			}
	}
}

MatrixPtr read_matrix(std::vector<std::string> lines)
{
	if(lines.empty())
		throw std::logic_error("Reading in matrix failed because lines vector is empty.");

	std::regex dynamic_regex("^dynamic ([0-9]+) ([0-9]+) ?\\[?([0-9,;\\.\\/]*)\\]?", std::regex_constants::icase);
	std::regex f_regex("^F_\\{?([0-9]+),([0-9]+)\\}?", std::regex_constants::icase);
	std::regex one_b_regex("^1_\\{?([0-9])\\}?", std::regex_constants::icase);
	std::regex terminator_regex("^T", std::regex_constants::icase);
	std::regex zero_regex("^0_\\{?([0-9]+),([0-9]+)\\}", std::regex_constants::icase);
	//std::regex diagonal_regex("^DIAG_\\{?([0-9]+),([0-1]+)\\}?", std::regex_constants::icase);
    std::regex diagonal_regex("^diagonal ([0-9]+) ?\\[?([0-9,\\.\\/]*)\\]?", std::regex_constants::icase);


    // dynamic matrix
	std::smatch matches;
	if(std::regex_match(lines[0], matches, dynamic_regex))
	{
		int n = std::stoi(matches[1]);
		int m = std::stoi(matches[2]);

		if(n > 63 || m > 63)
			throw std::logic_error("Dimensions of matrix too big.");

		std::string matrix_str = matches[3];

		MatrixPtr rtn = std::make_shared<DynamicMatrix>(n,m);
		auto& matrix = *rtn;
		// [1,2,3; 4,5,6; 7,8,9] format
		if(!matrix_str.empty())
		{
			boost::trim(matrix_str);

			std::vector<std::string> row_strs;
			boost::split(row_strs, matrix_str, boost::is_any_of(";"));
			if(row_strs.size() != (static_cast<std::size_t>(1) << m))
				throw std::logic_error(std::string("Rows: Provided array") + std::to_string(row_strs.size()) + " does not have right dimension " + std::to_string((1 << m)));

			unsigned long long i_row = 0;
			for(auto row_str : row_strs)
			{
				boost::trim_if(row_str,boost::is_any_of(" ,"));

				std::vector<std::string> value_strs;
				boost::split(value_strs, row_str, boost::is_any_of(","));

				if(value_strs.size() != (static_cast<std::size_t>(1) << n))
					throw std::logic_error(std::string("Cols: Provided array") + std::to_string(value_strs.size()) + " does not have right dimension " + std::to_string((1 << n)));

				unsigned long long i_col = 0;
				for(auto val : value_strs)
					matrix.set(i_row, i_col++, read_double(val));

				i_row++;
			}

		}

		if(!is_stochastic(*rtn))
			rtn->is_stochastic = false;

		return rtn;
	}

	// F matrix
	if(std::regex_match(lines[0], matches, f_regex))
	{
		int k = std::stoi(matches[1].str());
		int b = std::stoi(matches[2].str());

		return std::make_shared<FMatrix>(k,b);
	}

	// One b matrix
	if(std::regex_match(lines[0], matches, one_b_regex))
	{
		int b = std::stoi(matches[1].str());

		return std::make_shared<OneBMatrix>(b);
	}

	// terminator matrix
	if(std::regex_match(lines[0], matches, terminator_regex))
	{
		return std::make_shared<TerminatorMatrix>();
	}

	// F matrix
	if(std::regex_match(lines[0], matches, zero_regex))
	{
		int n = std::stoi(matches[1].str());
		int m = std::stoi(matches[2].str());

		return std::make_shared<ZeroMatrix>(n,m);
	}

    // diagonal matrix
    if(std::regex_match(lines[0], matches, diagonal_regex))
    {
        int k = std::stoi(matches[1]);

        if(k > 63)
            throw std::logic_error("Dimensions of matrix too big.");

        std::string diag_str = matches[2];

        MatrixPtr rtn = std::make_shared<DiagonalMatrix>(k);
        auto& matrix = *rtn;
        // [1,2,3,4] format
        if(!diag_str.empty())
        {
            boost::trim(diag_str);

            std::vector<std::string> val_strs;
            boost::split(val_strs, diag_str, boost::is_any_of(","));
            if(val_strs.size() != (static_cast<std::size_t>(1) << k))
                throw std::logic_error(std::string("Provided array") + std::to_string(val_strs.size()) + " does not have right dimension " + std::to_string((1 << k)));

            unsigned long long i_col = 0;
            for(auto val : val_strs){
                boost::trim_if(val,boost::is_any_of(" ,"));
                matrix.set(i_col, i_col, read_double(val));
                i_col++;
            }
        }

        if(!is_stochastic(*rtn))
            rtn->is_stochastic = false;

        return rtn;
    }

	throw std::logic_error(std::string("Matrix could not be read. Not able to parse line")+ "'" + lines[0] + "'.");
}

void print_matrix(std::ostream& ostr, const Matrix& matrix)
{
	ostr << "n: " << matrix.n << " m: " << matrix.m << std::endl;
	unsigned long long i_max_row = 1;
	unsigned long long i_max_col = 1;
	i_max_col = i_max_col << matrix.n;
	i_max_row = i_max_row << matrix.m;

	ostr.precision(4);


	ostr << std::setw(matrix.m+3) << " ";
	for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
		ostr << std::setw(6) << BitVec(i_col).to_string().substr(MAX_PLACES-matrix.n) << " ";	
	ostr << std::endl;

	ostr << std::setw(matrix.m+3) << " ";
	for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
		ostr << std::setw(5) << std::fixed << "------" << " ";	
	ostr << std::endl;

	for(unsigned long long i_row = 0; i_row < i_max_row; i_row++)
	{
		ostr << BitVec(i_row).to_string().substr(MAX_PLACES-matrix.m) << " | ";
		for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
		{
			ostr << std::setw(5) << std::fixed << matrix.get(i_row, i_col) << " ";	
		}
		ostr << "\n";	
	}

	ostr << std::endl;
}
