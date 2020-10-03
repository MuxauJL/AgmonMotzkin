#include "AgmonMotzkin.h"
#include <Debug.h>
#include <numeric>

bool AgmonMotzkin::init(const ClpSimplex& rModel)
{
	const double* pObjCoefficients = rModel.getObjCoefficients();
	const double* pRowLowers = rModel.getRowLower();
	const double* pRowUppers = rModel.getRowUpper();
	const double* pColLowers = rModel.getColLower();
	const double* pColUppers = rModel.getColUpper();
	int iNumRows = rModel.getNumRows();
	int iNumCols = rModel.getNumCols();

	optimizationDirection_ = rModel.optimizationDirection();
	objCoefficients_.assign(pObjCoefficients, pObjCoefficients + iNumCols);
	rowLowers_.assign(pRowLowers, pRowLowers + iNumRows);
	colLowers_.assign(pColLowers, pColLowers + iNumCols);
	auto* pMatrix = rModel.matrix();
	if (pMatrix)
	{
		// we store the row ordered matrix_
		if (pMatrix->isColOrdered())
		{
			matrix_.reverseOrderedCopyOf(*pMatrix);
		}
		else
		{
			matrix_ = *pMatrix;
		}
		MV_ASSERT(matrix_.isColOrdered() == false);
	}
	else
	{
		MV_LOGD("The model must contain a matrix");
		return false;
	}

	x_.resize(iNumCols, 0);
	coefficientSquaresSums_.clear();

	std::vector<double> rowUppers(pRowUppers, pRowUppers + iNumRows);
	convertToStandardFormMin(rowUppers);

	return true;
}

bool AgmonMotzkin::nextApproximation(double dEpsilon)
{
	MV_ASSERT(matrix_.isColOrdered() == false);
	MV_ASSERT(x_.size() == matrix_.getNumCols());

	auto [iRow, dStep] = prepareStep();
	if ((iRow < 0) || (dStep < dEpsilon))
	{
		return false;
	}
	MV_ASSERT(iRow >= 0);
	MV_ASSERT(dStep > 0);

	int iNumRows = matrix_.getNumRows();
	int iNumCols = matrix_.getNumCols();
	MV_ASSERT(iRow < iNumRows + iNumCols);
	// Indices iRow = iNumRows, ..., iNumRows + iNumCols - 1
	// correspond to constraints: 
	// x[iCol] >= 0, iCol = 0, ..., iNumCols - 1
	if (iRow >= iNumRows)
	{
		x_[iRow - iNumRows] += dStep;
		return true;
	}

	auto direction = matrix_.getVector(iRow);
	const double* pCoefficients = direction.getElements();
	const int iLength = direction.getNumElements();
	const int* pIndices = direction.getIndices();
	for (int iIndex = 0; iIndex < iLength; ++iIndex)
	{
		x_[pIndices[iIndex]] += dStep * pCoefficients[iIndex];
	}
	return true;
}

void AgmonMotzkin::convertToStandardFormMin(std::vector<double>& rRowUppers)
{
	double* pCoefficients = matrix_.getMutableElements();
	int* pStarts = matrix_.getMutableVectorStarts();
	int* pLengths = matrix_.getMutableVectorLengths();
	int* pIndices = matrix_.getMutableIndices();
	int iNumRows = matrix_.getNumRows();
	for (int iRow = 0; iRow < iNumRows; ++iRow)
	{
		// skip constraints: { lower <= AX <= +inf }
		if (rRowUppers[iRow] == COIN_DBL_MAX)
		{
			continue;
		}
		else
		{
			// convert constraints: { -inf <= AX <= upper } -> { -AX >= -upper }
			if (rowLowers_[iRow] == -COIN_DBL_MAX)
			{
				int iLength = pLengths[iRow];
				int iStart = pStarts[iRow];
				for (int iIndex = 0; iIndex < iLength; ++iIndex)
				{
					pCoefficients[iStart + iIndex] *= -1;
				}
				rowLowers_[iRow] = -rRowUppers[iRow];
				rRowUppers[iRow] = COIN_DBL_MAX;
			}
			else // convert constraints: { lower <= AX <= upper } -> { AX >= lower; -AX >= -upper }
			{
				int iLength = pLengths[iRow];
				MV_ASSERT(iLength >= 0);
				int iStart = pStarts[iRow];
				double* pCoefficientsStart = &pCoefficients[iStart];
				std::vector<double> negativeCoefficients;
				std::transform(pCoefficientsStart, pCoefficientsStart + iLength,
					std::back_inserter(negativeCoefficients),
					[](const double& elem)
					{
						return -elem;
					});
				MV_ASSERT(iLength == negativeCoefficients.size());				
				appendConstraint(iLength, &pIndices[iStart], negativeCoefficients.data(),
					-rRowUppers[iRow], COIN_DBL_MAX);

				rRowUppers[iRow] = COIN_DBL_MAX;

				pCoefficients = matrix_.getMutableElements();
				pStarts = matrix_.getMutableVectorStarts();
				pLengths = matrix_.getMutableVectorLengths();
				pIndices = matrix_.getMutableIndices();
			}
		}
	}
}

double AgmonMotzkin::coefficientSquaresSum(int iRowIndex)
{
	MV_ASSERT(matrix_.isColOrdered() == false);
	MV_ASSERT(iRowIndex >= 0);
	MV_ASSERT(iRowIndex < matrix_.getNumRows());

	auto it = coefficientSquaresSums_.find(iRowIndex);
	if (it != coefficientSquaresSums_.end())
	{
		return it->second;
	}

	double dResult = 0;
	auto coefficients = matrix_.getVector(iRowIndex);
	const double* pCoefficients = coefficients.getElements();
	const int iLength = coefficients.getNumElements();
	const int* pIndices = coefficients.getIndices();
	for (int iIndex = 0; iIndex < iLength; ++iIndex)
	{
		double dCoeff = coefficients[pIndices[iIndex]];
		dResult += dCoeff * dCoeff;
	}
	coefficientSquaresSums_[iRowIndex] = dResult;
	return dResult;
}

std::pair<int, double> AgmonMotzkin::prepareStep()
{
	MV_ASSERT(matrix_.isColOrdered() == false);

	double dMaxViolation = 0;
	int iIdx = -1;

	int iNumRows = matrix_.getNumRows();
	for (int iRow = 0; iRow < iNumRows; ++iRow)
	{
		double dResult = 0;
		auto coefficients = matrix_.getVector(iRow);
		const double* pCoefficients = coefficients.getElements();
		const int iLength = coefficients.getNumElements();
		const int* pIndices = coefficients.getIndices();
		for (int iIndex = 0; iIndex < iLength; ++iIndex)
		{
			dResult += pCoefficients[iIndex] * x_[pIndices[iIndex]];
		}
		double dLowerBound = rowLowers_[iRow];
		if (dResult < dLowerBound)
		{
			double dStep = (dLowerBound - dResult) / (coefficientSquaresSum(iRow));
			if (dMaxViolation < dStep)
			{
				dMaxViolation = dStep;
				iIdx = iRow;
			}
		}
	}

	// check that x[iCol] >= 0, iCol = 0, 1, ..., iNumCols - 1
	int iNumCols = matrix_.getNumCols();
	for (int iCol = 0; iCol < iNumCols; ++iCol)
	{
		double dStep = -x_[iCol];
		if (dMaxViolation < dStep)
		{
			dMaxViolation = dStep;
			iIdx = iNumRows + iCol;
		}
	}
	
	return { iIdx, dMaxViolation };
}

bool AgmonMotzkin::setCurrentPoint(const std::vector<double>& rCurrentPoint)
{
	size_t nSize = x_.size();
	if (nSize != rCurrentPoint.size())
	{
		MV_LOGD("The point size must be the same as the current: %d", nSize);
		return false;
	}
	x_ = rCurrentPoint;
	return true;
}

const std::vector<double>& AgmonMotzkin::getCurrentPoint()
{
	return x_;
}

double AgmonMotzkin::getCurrentObjectiveValue()
{
	int iNumberColumns = matrix_.getNumCols();
	double dResult = 0;
	for (int iCol = 0; iCol < iNumberColumns; ++iCol)
	{
		dResult += objCoefficients_[iCol] * x_[iCol];
	}
	return dResult;
}

void AgmonMotzkin::appendConstraint(int iNumCoefficients,
	const int* pIndices, const double* pCoefficients, 
	double dLowerBound, double dUpperBound)
{
	CoinPackedVector row(iNumCoefficients, pIndices, pCoefficients, false);
	matrix_.appendRow(row);
	rowLowers_.emplace_back(dLowerBound);
	// convert constraints: { AX <= upper } -> { -AX >= -upper }
	if (dUpperBound < COIN_DBL_MAX)
	{
		std::vector<double> negativeCoefficients;
		std::transform(pCoefficients, pCoefficients + iNumCoefficients,
			std::back_inserter(negativeCoefficients),
			[](const double& elem)
			{
				return -elem;
			});
		MV_ASSERT(iNumCoefficients == negativeCoefficients.size());
		CoinPackedVector negativeRow(iNumCoefficients, pIndices, negativeCoefficients.data(), false);
		matrix_.appendRow(negativeRow);
		rowLowers_.emplace_back(-dUpperBound);
	}
}

void AgmonMotzkin::printMatrix()
{
	MV_ASSERT(matrix_.isColOrdered() == false);

	int iNumberRows = matrix_.getNumRows();
	int iNumberCols = matrix_.getNumCols();
	for (int iRow = 0; iRow < iNumberRows; ++iRow)
	{
		printf("%d row:\n%f <= ", iRow, rowLowers_[iRow]);
		const auto row = matrix_.getVector(iRow);
		const double* pCoefficients = row.getElements();
		const int iLength = row.getNumElements();
		const int* pIndices = row.getIndices();
		for (int iIndex = 0; iIndex < iLength; ++iIndex)
		{
			printf("%f ", row[pIndices[iIndex]]);
		}
		printf("<= +inf\n");
	}
}
