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


	PINVOKE void Inflate(
		const Point3d* inputPoints,
		int inputCount,
		double  delta,
		Point3d** outputPoints0,
		int& outputCount0,
		int** offsetOriginalIndices,
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt
	)
	{
		Path64 outPath;
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);

		float miter_limit = 1.0;


		Paths64 solution = Clipper2Lib::InflatePaths(
			Paths64{ inputPath }, // expects Paths64
			delta * 1e6,        // scale delta as in other usages
			joinType,
			endType,
			miter_limit,
			arc_tolerance// assuming this is miter limit or arc tolerance
		);

		// Convert the first output path to Point3d array
		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount0 = static_cast<int>(solution[0].size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(solution[0], *outputPoints0, outputCount0);

			std::vector<int> offsetOriginalIndicesVec;
			*offsetOriginalIndices = new int[outputCount0];
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
		}
	}
	PINVOKE void InflateVariable(
		const Point3d* inputPoints,
		int inputCount,
		double* deltas,
		Point3d** outputPoints0,
		int& outputCount0,
		int** offsetOriginalIndices, 
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt
	)
	{
		Path64 outPath;
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);


		float miter_limit = 1.0;
		ClipperOffset co(miter_limit, arc_tolerance);

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
		co.Execute(1e6, solution);
		// Convert the first output path to Point3d array
		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount0 = static_cast<int>(solution[0].size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(solution[0], *outputPoints0, outputCount0);

			std::vector<int> offsetOriginalIndicesVec;
			*offsetOriginalIndices = new int[outputCount0];
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
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt
	)
	{
		Path64 outPath;
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);
		Paths64 subject;



		subject.push_back(inputPath);

		float miter_limit = 1.0;
		ClipperOffset co(miter_limit, arc_tolerance);

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
		
		co.Execute(1e6, solution);
		SimplifyPaths(solution, 1 * 1e6);

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