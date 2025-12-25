#pragma once
#define PINVOKE extern "C" __declspec(dllexport)

#include "clipper2/clipper.offset.h"
#include "clipper2/clipper.h"

namespace Clipper2Lib {

	const double floating_point_tolerance = 1e-12;
	const double default_arc_tolerance = 50;
	const double M_PI = 3.14159265358979323846;

	struct Point3d {
		double x;
		double y;
		double z;

		Point3d() : x(0), y(0), z(0) {}
		Point3d(double x_, double y_, double z_)
			: x(x_), y(y_), z(z_) {
		}
	};	Path64 ConvertToPath64(const Point3d* points, int count) {
		Path64 path;
		path.reserve(count);
		for (int i = 0; i < count; i++) {
			// Clipper2 uses integer coordinates, so we scale the doubles
			// Using 1e6 as scale factor for 6 decimal places of precision
			const double scale = 1e6;
			int64_t x = static_cast<int64_t>(points[i].x * scale);
			int64_t y = static_cast<int64_t>(points[i].y * scale);
			int64_t z = static_cast<int64_t>(points[i].z * scale);
			path.push_back(Point64(x, y, z));
		}
		return path;
	}
	static inline Point64 GetPerpendic(const Point64& pt, const PointD& norm, double delta)
	{
		return Point64(pt.x + norm.x * delta, pt.y + norm.y * delta, pt.z, pt.w, pt.o);
	}

	void ConvertFromPath64(const Path64& path, Point3d* points, int& count) {
		const double scale = 1e-6; // Reverse the scaling
		count = 0;
		for (const Point64& p : path) {
			points[count].x = p.x * scale;
			points[count].y = p.y * scale;
			points[count].z = p.z * scale;
			count++;
		}
	}

	PINVOKE void DeletePoints(Point3d* points) {
		delete[] points;
	}

	void SmoothZValuesAcrossPlateaus(Path64& path) {
		if (path.size() < 3) return;

		std::vector<int64_t> originalZ(path.size());
		for (size_t i = 0; i < path.size(); i++) {
			originalZ[i] = path[i].z;
		}

		// Find plateau boundaries (first occurrence of each unique z value)
		std::vector<std::pair<size_t, int64_t>> keyPoints; // (index, z-value)
		keyPoints.push_back({ 0, originalZ[0] });

		for (size_t i = 1; i < path.size(); i++) {
			if (originalZ[i] != originalZ[i - 1]) {
				keyPoints.push_back({ i, originalZ[i] });
			}
		}

		// Add the last point if it's not already a key point
		if (keyPoints.back().first != path.size() - 1) {
			keyPoints.push_back({ path.size() - 1, originalZ.back() });
		}

		// Interpolate across the entire path using key points
		for (size_t i = 0; i < path.size(); i++) {
			// Find which segment this index falls into
			for (size_t j = 0; j < keyPoints.size() - 1; j++) {
				size_t startIdx = keyPoints[j].first;
				size_t endIdx = keyPoints[j + 1].first;

				if (i >= startIdx && i <= endIdx) {
					int64_t startZ = keyPoints[j].second;
					int64_t endZ = keyPoints[j + 1].second;
					double t = static_cast<double>(i - startIdx) / static_cast<double>(endIdx - startIdx);
					path[i].z = static_cast<int64_t>(startZ + t * (endZ - startZ));
					break;
				}
			}
		}
	}
	bool IsPathClosed(const Point3d* points, int count) {
		if (count < 2) return false;

		const Point3d& first = points[0];
		const Point3d& last = points[count - 1];

		// Compare with small epsilon for floating point tolerance
		const double epsilon = 1e-9;
		return (std::abs(first.x - last.x) < epsilon &&
			std::abs(first.y - last.y) < epsilon &&
			std::abs(first.z - last.z) < epsilon);
	}

	PINVOKE void Inflate(
		const Point3d* inputPoints,
		int inputCount,
		double  delta,
		Point3d** outputPoints0,
		int& outputCount0,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt,
		double epsilonForSimplifying = 10.0,
		bool simplifyBeforeOffset = true,
		bool smoothZ = true,
		bool simplifyPath = true
	)
	{
		Path64 outPath;
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);
		bool isClosedPath = IsPathClosed(inputPoints, inputCount);
		if (isClosedPath)
			if (endType != EndType::Round)
				endType = EndType::Polygon;

		if (simplifyBeforeOffset) {
			inputPath = Clipper2Lib::SimplifyPath(inputPath, epsilonForSimplifying, isClosedPath);
		}

		Paths64 outPaths = Clipper2Lib::InflatePaths(
			Paths64{ inputPath }, // expects Paths64
			delta * 1e6,        // scale delta as in other usages
			joinType,
			endType,
			3,
			100// assuming this is miter limit or arc tolerance
		);

		// Convert the first output path to Point3d array
		if (!outPaths.empty()) {
			outPaths[0].push_back(outPaths[0].front());
			outputCount0 = static_cast<int>(outPaths[0].size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(outPaths[0], *outputPoints0, outputCount0);
		}
		else {
			outputCount0 = 0;
			*outputPoints0 = nullptr;
		}
	}
	PINVOKE void InflateVariable(
		const Point3d* inputPoints,
		int inputCount,
		double* deltas,
		Point3d** outputPoints0,
		int& outputCount0,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt,
		double epsilonForSimplifying = 10.0,
		bool simplifyBeforeOffset = true,
		bool smoothZ = true,
		bool simplifyPath = true
	)
	{
		Path64 outPath;
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);
		bool isClosedPath = IsPathClosed(inputPoints, inputCount);

		if (simplifyBeforeOffset) {
			inputPath = Clipper2Lib::SimplifyPath(inputPath, epsilonForSimplifying, isClosedPath);
		}

		ClipperOffset co;
		Paths64 subject;
		subject.push_back(inputPath);

		co.AddPaths(subject, joinType, endType);

		co.SetDeltaCallback([deltas](const Path64& path,
			const PathD& path_norms, size_t curr_idx, size_t prev_idx)
			{
				return  deltas[curr_idx] * 1e6;
			});

		//  solution
		Paths64 solution;
		co.Execute(1.0, solution);

		// Convert the first output path to Point3d array
		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount0 = static_cast<int>(solution[0].size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(solution[0], *outputPoints0, outputCount0);
		}
		else {
			outputCount0 = 0;
			*outputPoints0 = nullptr;
		}
	}
	struct DeltaSelector
	{
		const bool* isAnotherSide;
		const bool* isEnding;        // ?? NEW
		double* left;
		double* right;

		double operator()(const Path64&,
			const PathD&,
			size_t curr_idx,
			size_t) const
		{
			if (*isEnding)
			{
				// swapped behavior
				return (*isAnotherSide ?
					right[curr_idx] :
					left[curr_idx]) * 1e6;
			}
			else
			{
				// default behavior
				return (*isAnotherSide ?
					left[curr_idx] :
					right[curr_idx]) * 1e6;
			}
		}
	};

	PINVOKE void InflateVariableBothSide(
		const Point3d* inputPoints,
		int inputCount,
		double* leftSideDeltas,
		double* rightSideDeltas,
		Point3d** outputPoints0,
		int& outputCount0,
		int** wIndices,
		int** offsetOriginalIndices,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt,
		double epsilonForSimplifying = 10.0,
		bool simplifyBeforeOffset = true,
		bool smoothZ = true,
		bool simplifyPath = true
	)
	{
		Path64 outPath;
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);
		bool isClosedPath = IsPathClosed(inputPoints, inputCount);
		if (simplifyBeforeOffset) {
			inputPath = Clipper2Lib::SimplifyPath(inputPath, epsilonForSimplifying, isClosedPath);
		}
		ClipperOffset co;
		Paths64 subject;
		subject.push_back(inputPath);
		co.AddPaths(subject, joinType, endType);
		DeltaSelector selector{
			&co.is_another_side,
			&co.ending_flag,
			leftSideDeltas,
			rightSideDeltas
		};

		co.SetDeltaCallback(selector);
		// solution
		Paths64 solution;
		co.Execute(1.0, solution);

		// Convert the first output path to Point3d array
		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount0 = static_cast<int>(solution[0].size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(solution[0], *outputPoints0, outputCount0);

			// Collect all  w values
			std::vector<int> wTrueIndicesVec;
			std::vector<int> offsetOriginalIndicesVec;
			*wIndices = new int[outputCount0];
			*offsetOriginalIndices = new int[outputCount0];

			for (int i = 0; i < static_cast<int>(solution[0].size()); ++i)
			{
				wTrueIndicesVec.push_back(solution[0][i].w);

			}

			// Allocate and copy indices
			std::copy(wTrueIndicesVec.begin(), wTrueIndicesVec.end(), *wIndices);
			// Collect all offset original indices

			for (int i = 0; i < static_cast<int>(solution[0].size()); ++i) {
				offsetOriginalIndicesVec.push_back(solution[0][i].o);
			}

			// Allocate and copy to output pointer
			std::copy(offsetOriginalIndicesVec.begin(),
				offsetOriginalIndicesVec.end(),
				*offsetOriginalIndices);

		}
		else {
			outputCount0 = 0;
			*outputPoints0 = nullptr;
			*wIndices = nullptr;
			*offsetOriginalIndices = nullptr;
		}
	}
	extern "C" __declspec(dllexport)
		void FreeMemory(Point3d* ptr)
	{
		delete[] ptr;
	}


#include <vector>
#include <cmath>




}