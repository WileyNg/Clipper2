#pragma once
#define PINVOKE extern "C" __declspec(dllexport)

#include "clipper2/clipper.offset.h"
#include "clipper2/clipper.h"
#include <combaseapi.h> 

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
	};	
	Path64 ConvertToPath64(const Point3d* points, int count) {
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
	Path64 ConvertToPath64WithPathIndex(const Point3d* points, int count, int pathIndex, const int* wList, const int* oList) {
		Path64 path;
		path.reserve(count);
		for (int i = 0; i < count; i++) {
			const double scale = 1e6;
			int64_t x = static_cast<int64_t>(points[i].x * scale);
			int64_t y = static_cast<int64_t>(points[i].y * scale);
			int64_t z = static_cast<int64_t>(points[i].z * scale);

			int w = (wList != nullptr) ? wList[i] : 0;
			int o = (oList != nullptr) ? oList[i] : 0;

			Point64 pt;
			pt.Init(x, y, z, w, o, pathIndex);
			path.push_back(pt);
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
	// Standalone Z callback function for Union operation
	static void UnionZCB(const Point64& bot1, const Point64& top1,
		const Point64& bot2, const Point64& top2, Point64& ip)
	{
		// Follow the same logic as ClipperOffset::ZCB
		if (bot1.z && ((bot1.z == bot2.z) || (bot1.z == top2.z)))
			ip.z = bot1.z;
		else if (bot2.z && (bot2.z == top1.z))
			ip.z = bot2.z;
		else if (top1.z && (top1.z == top2.z))
			ip.z = top1.z;
		// No external zCallback64_ in standalone version, so default to top1.z
		else
			ip.z = top1.z;

		// FIXED: Only set w, o, p_i if ip doesn't already have valid metadata
		if (ip.p_i == -1 || ip.o == -1 || ip.w == -1)
		{
			// ip needs metadata - choose from a valid source (prefer top1)
			if (top1.p_i != -1)
			{
				ip.w = top1.w;
				ip.o = top1.o;
				ip.p_i = top1.p_i;
			}
			else if (bot1.p_i != -1)
			{
				ip.w = bot1.w;
				ip.o = bot1.o;
				ip.p_i = bot1.p_i;
			}
			else if (top2.p_i != -1)
			{
				ip.w = top2.w;
				ip.o = top2.o;
				ip.p_i = top2.p_i;
			}
			else if (bot2.p_i != -1)
			{
				ip.w = bot2.w;
				ip.o = bot2.o;
				ip.p_i = bot2.p_i;
			}
		}
	}

	PINVOKE void Union(
		const Point3d** inputPolylines,
		const int* inputCounts,
		int polylineCount,
		int** inputWList,
		int** inputOList,
		Point3d*** outputPolylines,
		int** outputCounts,
		int& outputPolylineCount,
		int*** wIndices,
		int*** offsetOriginalIndices,
		int*** pathIndices,
		bool preserveCollinear = false,
		FillRule fillRule = FillRule::Positive
	)
	{
		// Convert all input polylines to Path64
		Paths64 subjects;
		subjects.reserve(polylineCount);
		for (int i = 0; i < polylineCount; ++i)
		{
			const int* wList = (inputWList != nullptr) ? inputWList[i] : nullptr;
			const int* oList = (inputOList != nullptr) ? inputOList[i] : nullptr;
			Path64 path = ConvertToPath64WithPathIndex(inputPolylines[i], inputCounts[i], i, wList, oList);
			subjects.push_back(path);
		}

		// Perform union operation
		Paths64 solution;
		Clipper64 clipper;
		clipper.PreserveCollinear(preserveCollinear);
#ifdef USINGZ
		clipper.SetZCallback(UnionZCB);
#endif
		clipper.AddSubject(subjects);
		clipper.Execute(ClipType::Union, fillRule, solution);

		// Convert results back to Point3d arrays
		if (!solution.empty())
		{
			outputPolylineCount = static_cast<int>(solution.size());

			// Allocate outer arrays using CoTaskMemAlloc
			*outputPolylines = (Point3d**)CoTaskMemAlloc(sizeof(Point3d*) * outputPolylineCount);
			*outputCounts = (int*)CoTaskMemAlloc(sizeof(int) * outputPolylineCount);
			*wIndices = (int**)CoTaskMemAlloc(sizeof(int*) * outputPolylineCount);
			*offsetOriginalIndices = (int**)CoTaskMemAlloc(sizeof(int*) * outputPolylineCount);
			*pathIndices = (int**)CoTaskMemAlloc(sizeof(int*) * outputPolylineCount);

			// Allocate and fill inner arrays
			for (int i = 0; i < outputPolylineCount; ++i)
			{
				int count = static_cast<int>(solution[i].size());
				(*outputCounts)[i] = count;

				// Allocate inner arrays using CoTaskMemAlloc
				(*outputPolylines)[i] = (Point3d*)CoTaskMemAlloc(sizeof(Point3d) * count);
				(*wIndices)[i] = (int*)CoTaskMemAlloc(sizeof(int) * count);
				(*offsetOriginalIndices)[i] = (int*)CoTaskMemAlloc(sizeof(int) * count);
				(*pathIndices)[i] = (int*)CoTaskMemAlloc(sizeof(int) * count);

				ConvertFromPath64(solution[i], (*outputPolylines)[i], count);

				// Extract w, o, and p_i indices
				for (int j = 0; j < count; ++j)
				{
					(*wIndices)[i][j] = solution[i][j].w;
					(*offsetOriginalIndices)[i][j] = solution[i][j].o;
					(*pathIndices)[i][j] = solution[i][j].p_i;
				}
			}
		}
		else
		{
			outputPolylineCount = 0;
			*outputPolylines = nullptr;
			*outputCounts = nullptr;
			*wIndices = nullptr;
			*offsetOriginalIndices = nullptr;
			*pathIndices = nullptr;
		}
	}
	PINVOKE void UnionAndIntersection(
		const Point3d** subjectPolylines,
		const int* subjectCounts,
		int subjectPolylineCount,
		int** subjectWList,              // NEW: Array of w indices for each subject polyline
		int** subjectOList,              // NEW: Array of o indices for each subject polyline
		const Point3d** clipPolylines,
		const int* clipCounts,
		int clipPolylineCount,
		int** clipWList,                 // NEW: Array of w indices for each clip polyline
		int** clipOList,                 // NEW: Array of o indices for each clip polyline
		Point3d*** outputPolylines,      // Array of output polyline pointers
		int** outputCounts,              // Array of vertex counts for each output polyline
		int& outputPolylineCount,        // Number of output polylines
		int*** wIndices,                 // Array of w indices for each output polyline
		int*** offsetOriginalIndices,    // Array of o indices for each output polyline
		int*** pathIndices,              // Array of p_i (path index) for each output polyline
		bool preserveCollinear = false,
		FillRule fillRule = FillRule::Positive
	)
	{
		// Convert subject polylines to Path64
	// Convert subject polylines to Path64
		Paths64 subjects;
		subjects.reserve(subjectPolylineCount);
		for (int i = 0; i < subjectPolylineCount; ++i)
		{
			const int* wList = (subjectWList != nullptr) ? subjectWList[i] : nullptr;
			const int* oList = (subjectOList != nullptr) ? subjectOList[i] : nullptr;
			Path64 path = ConvertToPath64WithPathIndex(subjectPolylines[i], subjectCounts[i], i, wList, oList);
			subjects.push_back(path);
		}

		// Convert clip polylines to Path64
		Paths64 clips;
		for (int i = 0; i < clipPolylineCount; ++i)
		{
			// Create arrays filled with -1 for clip polylines
			int count = clipCounts[i];
			int* wList = new int[count];
			int* oList = new int[count];

			// Fill with -1
			for (int j = 0; j < count; ++j)
			{
				wList[j] = -1;
				oList[j] = -1;
			}

			Path64 path = ConvertToPath64WithPathIndex(clipPolylines[i], clipCounts[i], -1, wList, oList);
			clips.push_back(path);

			// Clean up temporary arrays
			delete[] wList;
			delete[] oList;
		}

		// Perform union operation first
		Paths64 unionSolution;
		Clipper64 clipper;
		clipper.PreserveCollinear(preserveCollinear);
#ifdef USINGZ
		clipper.SetZCallback(UnionZCB);
#endif
		clipper.AddSubject(subjects);
		clipper.Execute(ClipType::Union, fillRule, unionSolution);

		// Now perform intersection on the union result with the clip polylines
		Paths64 solution;
		clipper.Clear();
		clipper.PreserveCollinear(preserveCollinear);
#ifdef USINGZ
		clipper.SetZCallback(UnionZCB);
#endif
		clipper.AddSubject(unionSolution);
		clipper.AddClip(clips);
		clipper.Execute(ClipType::Intersection, fillRule, solution);

		// Convert results back to Point3d arrays
		if (!solution.empty())
		{
			outputPolylineCount = static_cast<int>(solution.size());

			// Allocate outer arrays using CoTaskMemAlloc
			*outputPolylines = (Point3d**)CoTaskMemAlloc(sizeof(Point3d*) * outputPolylineCount);
			*outputCounts = (int*)CoTaskMemAlloc(sizeof(int) * outputPolylineCount);
			*wIndices = (int**)CoTaskMemAlloc(sizeof(int*) * outputPolylineCount);
			*offsetOriginalIndices = (int**)CoTaskMemAlloc(sizeof(int*) * outputPolylineCount);
			*pathIndices = (int**)CoTaskMemAlloc(sizeof(int*) * outputPolylineCount);

			// Allocate and fill inner arrays
			for (int i = 0; i < outputPolylineCount; ++i)
			{
				int count = static_cast<int>(solution[i].size());
				(*outputCounts)[i] = count;

				// Allocate inner arrays using CoTaskMemAlloc
				(*outputPolylines)[i] = (Point3d*)CoTaskMemAlloc(sizeof(Point3d) * count);
				(*wIndices)[i] = (int*)CoTaskMemAlloc(sizeof(int) * count);
				(*offsetOriginalIndices)[i] = (int*)CoTaskMemAlloc(sizeof(int) * count);
				(*pathIndices)[i] = (int*)CoTaskMemAlloc(sizeof(int) * count);

				ConvertFromPath64(solution[i], (*outputPolylines)[i], count);

				// Extract w, o, and p_i indices
				for (int j = 0; j < count; ++j)
				{
					(*wIndices)[i][j] = solution[i][j].w;
					(*offsetOriginalIndices)[i][j] = solution[i][j].o;
					(*pathIndices)[i][j] = solution[i][j].p_i;
				}
			}
		}
		else
		{
			outputPolylineCount = 0;
			*outputPolylines = nullptr;
			*outputCounts = nullptr;
			*wIndices = nullptr;
			*offsetOriginalIndices = nullptr;
			*pathIndices = nullptr;
		}
	}

	PINVOKE void InflateVariableBothSideWithoutUnion(
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

		co.ExecuteWithoutUnion(1e6, solution);

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