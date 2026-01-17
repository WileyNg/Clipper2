#pragma once
#define PINVOKE extern "C" __declspec(dllexport)

#include "clipper2/clipper.offset.h"
#include "clipper2/clipper.h"
#include <combaseapi.h> 

namespace Clipper2Lib {

	const double floating_point_tolerance = 1e-12;
	const double default_arc_tolerance = 50;
	const double M_PI = 3.14159265358979323846;

	// Simple Point3d for geometry only
	struct Point3d {
		double x;
		double y;
		double z;

		Point3d() : x(0), y(0), z(0) {}
		Point3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
	};

	// NEW: BoundaryPoint bundles geometry + metadata
	struct BoundaryPoint {
		Point3d point;
		int side;          // w - Position enum as int (0=Ending, 1=Left, 2=Right)
		int originalIndex; // o - index in original polyline
		int pathIndex;     // p_i - which path this came from

		BoundaryPoint() : point(), side(-1), originalIndex(-1), pathIndex(-1) {}
		BoundaryPoint(Point3d p, int s = -1, int o = -1, int pi = -1)
			: point(p), side(s), originalIndex(o), pathIndex(pi) {
		}
	};

	// ============================================================================
	// CONVERSION HELPERS
	// ============================================================================

	Path64 BoundaryPointsToPath64(const BoundaryPoint* points, int count) {
		Path64 path;
		path.reserve(count);
		const double scale = 1e6;

		for (int i = 0; i < count; i++) {
			Point64 pt;
			pt.Init(
				points[i].point.x * scale,
				points[i].point.y * scale,
				points[i].point.z * scale,
				points[i].side,
				points[i].originalIndex,
				points[i].pathIndex
			);
			path.push_back(pt);
		}
		return path;
	}

	void Path64ToBoundaryPoints(const Path64& path, BoundaryPoint* points, int& count) {
		const double scale = 1e-6;
		count = 0;
		for (const Point64& p : path) {
			points[count].point.x = p.x * scale;
			points[count].point.y = p.y * scale;
			points[count].point.z = p.z * scale;
			points[count].side = p.w;
			points[count].originalIndex = p.o;
			points[count].pathIndex = p.p_i;
			count++;
		}
	}

	// Legacy conversions (kept for backward compatibility)
	Path64 ConvertToPath64(const Point3d* points, int count) {
		Path64 path;
		path.reserve(count);
		const double scale = 1e6;
		for (int i = 0; i < count; i++) {
			int64_t x = static_cast<int64_t>(points[i].x * scale);
			int64_t y = static_cast<int64_t>(points[i].y * scale);
			int64_t z = static_cast<int64_t>(points[i].z * scale);
			path.push_back(Point64(x, y, z));
		}
		return path;
	}

	Path64 ConvertToPath64WithData(const Point3d* points, int count, const int* pathIndex, const int* oList) {
		Path64 path;
		path.reserve(count);
		const double scale = 1e6;
		for (int i = 0; i < count; i++) {
			int64_t x = static_cast<int64_t>(points[i].x * scale);
			int64_t y = static_cast<int64_t>(points[i].y * scale);
			int64_t z = static_cast<int64_t>(points[i].z * scale);

			int w = -1;
			int o = (oList != nullptr) ? oList[i] : 0;
			int p_i = (pathIndex != nullptr) ? pathIndex[i] : -1;
			Point64 pt;
			pt.Init(x, y, z, w, o, p_i);
			path.push_back(pt);
		}
		return path;
	}

	void ConvertFromPath64(const Path64& path, Point3d* points, int& count) {
		const double scale = 1e-6;
		count = 0;
		for (const Point64& p : path) {
			points[count].x = p.x * scale;
			points[count].y = p.y * scale;
			points[count].z = p.z * scale;
			count++;
		}
	}

	// ============================================================================
	// MEMORY MANAGEMENT
	// ============================================================================

	PINVOKE void FreeBoundaryPoints(BoundaryPoint* ptr) {
		CoTaskMemFree(ptr);
	}

	PINVOKE void FreeBoundaryPointsArray(BoundaryPoint** ptr, int count) {
		if (ptr) {
			for (int i = 0; i < count; i++) {
				if (ptr[i]) CoTaskMemFree(ptr[i]);
			}
			CoTaskMemFree(ptr);
		}
	}

	PINVOKE void DeletePoints(Point3d* points) {
		delete[] points;
	}

	PINVOKE void FreeMemory(Point3d* ptr) {
		delete[] ptr;
	}

	// ============================================================================
	// DELTA SELECTOR (for both-side operations)
	// ============================================================================

	struct DeltaSelector {
		const bool* isAnotherSide;
		const bool* isEnding;
		double* left;
		double* right;

		double operator()(const Path64&, const PathD&, size_t curr_idx, size_t) const {
			if (*isEnding) {
				return (*isAnotherSide ? right[curr_idx] : left[curr_idx]) * 1e6;
			}
			else {
				return (*isAnotherSide ? left[curr_idx] : right[curr_idx]) * 1e6;
			}
		}
	};

	static void UnionZCB(const Point64& bot1, const Point64& top1,
		const Point64& bot2, const Point64& top2, Point64& ip)
	{
		if (bot1.z && ((bot1.z == bot2.z) || (bot1.z == top2.z)))
			ip.z = bot1.z;
		else if (bot2.z && (bot2.z == top1.z))
			ip.z = bot2.z;
		else if (top1.z && (top1.z == top2.z))
			ip.z = top1.z;
		else
			ip.z = top1.z;

		if (ip.p_i == -1 || ip.o == -1 || ip.w == -1) {
			if (top1.p_i != -1) {
				ip.w = top1.w;
				ip.o = top1.o;
				ip.p_i = top1.p_i;
			}
			else if (bot1.p_i != -1) {
				ip.w = bot1.w;
				ip.o = bot1.o;
				ip.p_i = bot1.p_i;
			}
			else if (top2.p_i != -1) {
				ip.w = top2.w;
				ip.o = top2.o;
				ip.p_i = top2.p_i;
			}
			else if (bot2.p_i != -1) {
				ip.w = bot2.w;
				ip.o = bot2.o;
				ip.p_i = bot2.p_i;
			}
		}
	}

	// ============================================================================
	// LEGACY FUNCTIONS (keep for backward compatibility if needed)
	// ============================================================================

	PINVOKE void Inflate(
		const Point3d* inputPoints,
		int inputCount,
		double delta,
		Point3d** outputPoints0,
		int& outputCount0,
		int** offsetOriginalIndices,
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt)
	{
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);
		float miter_limit = 1.0;

		Paths64 solution = Clipper2Lib::InflatePaths(
			Paths64{ inputPath },
			delta * 1e6,
			joinType,
			endType,
			miter_limit,
			arc_tolerance
		);

		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount0 = static_cast<int>(solution[0].size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(solution[0], *outputPoints0, outputCount0);

			*offsetOriginalIndices = new int[outputCount0];
			for (int i = 0; i < outputCount0; ++i) {
				(*offsetOriginalIndices)[i] = solution[0][i].o;
			}
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
		const int* pathIndices,
		Point3d** outputPoints0,
		int& outputCount0,
		int** outputOffsetOriginalIndices,
		int** outputPathIndices,
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt)
	{
		Path64 inputPath = ConvertToPath64WithData(inputPoints, inputCount, pathIndices, nullptr);
		float miter_limit = 1.0;
		ClipperOffset co(miter_limit, arc_tolerance);
		Paths64 subject;
		subject.push_back(inputPath);
		co.AddPaths(subject, joinType, endType);
		co.SetDeltaCallback([deltas](const Path64& path, const PathD& path_norms,
			size_t curr_idx, size_t prev_idx) {
				return deltas[curr_idx] * 1e6;
			});

		Paths64 solution;
		co.Execute(1e6, solution);

		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount0 = static_cast<int>(solution[0].size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(solution[0], *outputPoints0, outputCount0);

			*outputOffsetOriginalIndices = new int[outputCount0];
			*outputPathIndices = new int[outputCount0];
			for (int i = 0; i < outputCount0; ++i) {
				(*outputOffsetOriginalIndices)[i] = solution[0][i].o;
				(*outputPathIndices)[i] = solution[0][i].p_i;
			}
		}
		else {
			outputCount0 = 0;
			*outputPoints0 = nullptr;
			*outputOffsetOriginalIndices = nullptr;
			*outputPathIndices = nullptr;
		}
	}

	// ============================================================================
	// REFACTORED FUNCTIONS - Using BoundaryPoint (keep same names!)
	// ============================================================================

	PINVOKE void InflateVariableWithoutUnion(
		const BoundaryPoint* inputPoints,
		int inputCount,
		double* deltas,
		BoundaryPoint** outputPoints,
		int& outputCount,
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt)
	{
		Path64 inputPath = BoundaryPointsToPath64(inputPoints, inputCount);

		float miter_limit = 1.0;
		ClipperOffset co(miter_limit, arc_tolerance);
		Paths64 subject;
		subject.push_back(inputPath);
		co.AddPaths(subject, joinType, endType);

		co.SetDeltaCallback([deltas](const Path64& path, const PathD& path_norms,
			size_t curr_idx, size_t prev_idx) {
				return deltas[curr_idx] * 1e6;
			});

		Paths64 solution;
		co.ExecuteWithoutUnion(1e6, solution);

		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount = static_cast<int>(solution[0].size());

			*outputPoints = (BoundaryPoint*)CoTaskMemAlloc(sizeof(BoundaryPoint) * outputCount);
			Path64ToBoundaryPoints(solution[0], *outputPoints, outputCount);
		}
		else {
			outputCount = 0;
			*outputPoints = nullptr;
		}
	}

	PINVOKE void InflateVariableBothSide(
		const BoundaryPoint* inputPoints,
		int inputCount,
		double* leftSideDeltas,
		double* rightSideDeltas,
		BoundaryPoint** outputPoints,
		int& outputCount,
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt)
	{
		Path64 inputPath = BoundaryPointsToPath64(inputPoints, inputCount);
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

		Paths64 solution;
		co.Execute(1e6, solution);

		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount = static_cast<int>(solution[0].size());

			*outputPoints = (BoundaryPoint*)CoTaskMemAlloc(sizeof(BoundaryPoint) * outputCount);
			Path64ToBoundaryPoints(solution[0], *outputPoints, outputCount);
		}
		else {
			outputCount = 0;
			*outputPoints = nullptr;
		}
	}

	PINVOKE void InflateVariableBothSideWithoutUnion(
		const BoundaryPoint* inputPoints,
		int inputCount,
		double* leftSideDeltas,
		double* rightSideDeltas,
		BoundaryPoint** outputPoints,
		int& outputCount,
		double arc_tolerance,
		JoinType joinType = JoinType::Round,
		EndType endType = EndType::Butt)
	{
		Path64 inputPath = BoundaryPointsToPath64(inputPoints, inputCount);
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

		Paths64 solution;
		co.ExecuteWithoutUnion(1e6, solution);

		if (!solution.empty()) {
			solution[0].push_back(solution[0].front());
			outputCount = static_cast<int>(solution[0].size());

			*outputPoints = (BoundaryPoint*)CoTaskMemAlloc(sizeof(BoundaryPoint) * outputCount);
			Path64ToBoundaryPoints(solution[0], *outputPoints, outputCount);
		}
		else {
			outputCount = 0;
			*outputPoints = nullptr;
		}
	}

	PINVOKE void Union(
		const BoundaryPoint** inputPolylines,
		const int* inputCounts,
		int polylineCount,
		BoundaryPoint*** outputPolylines,
		int** outputCounts,
		int& outputPolylineCount,
		bool preserveCollinear = false,
		FillRule fillRule = FillRule::NonZero)
	{
		Paths64 subjects;
		subjects.reserve(polylineCount);

		for (int i = 0; i < polylineCount; ++i) {
			Path64 path = BoundaryPointsToPath64(inputPolylines[i], inputCounts[i]);
			if (!Clipper2Lib::IsPositive(path)) {
				std::reverse(path.begin(), path.end());
			}
			subjects.push_back(std::move(path));
		}

		Paths64 solution;
		Clipper64 clipper;
		clipper.PreserveCollinear(preserveCollinear);
#ifdef USINGZ
		clipper.SetZCallback(UnionZCB);
#endif
		clipper.AddSubject(subjects);
		clipper.Execute(ClipType::Union, fillRule, solution);

		if (!solution.empty()) {
			outputPolylineCount = static_cast<int>(solution.size());

			*outputPolylines = (BoundaryPoint**)CoTaskMemAlloc(sizeof(BoundaryPoint*) * outputPolylineCount);
			*outputCounts = (int*)CoTaskMemAlloc(sizeof(int) * outputPolylineCount);

			for (int i = 0; i < outputPolylineCount; ++i) {
				int count = static_cast<int>(solution[i].size());
				(*outputCounts)[i] = count;

				(*outputPolylines)[i] = (BoundaryPoint*)CoTaskMemAlloc(sizeof(BoundaryPoint) * count);
				Path64ToBoundaryPoints(solution[i], (*outputPolylines)[i], count);
			}
		}
		else {
			outputPolylineCount = 0;
			*outputPolylines = nullptr;
			*outputCounts = nullptr;
		}
	}

	PINVOKE void UnionAndIntersection(
		const BoundaryPoint** subjectPolylines,
		const int* subjectCounts,
		int subjectPolylineCount,
		const Point3d** clipPolylines,
		const int* clipCounts,
		int clipPolylineCount,
		BoundaryPoint*** outputPolylines,
		int** outputCounts,
		int& outputPolylineCount,
		bool preserveCollinear = false,
		FillRule fillRule = FillRule::NonZero)
	{
		// Convert subject polylines (with metadata)
		Paths64 subjects;
		subjects.reserve(subjectPolylineCount);
		for (int i = 0; i < subjectPolylineCount; ++i) {
			Path64 path = BoundaryPointsToPath64(subjectPolylines[i], subjectCounts[i]);
			if (!Clipper2Lib::IsPositive(path)) {
				std::reverse(path.begin(), path.end());
			}

			subjects.push_back(std::move(path));
		}

		// Convert clip polylines (simple geometry only)
		Paths64 clips;
		clips.reserve(clipPolylineCount);
		for (int i = 0; i < clipPolylineCount; ++i) {
			Path64 path = ConvertToPath64(clipPolylines[i], clipCounts[i]);
			clips.push_back(std::move(path));
		}

		// Union subjects first
		Paths64 unionSolution;
		Clipper64 clipper;
		clipper.PreserveCollinear(preserveCollinear);
#ifdef USINGZ
		clipper.SetZCallback(UnionZCB);
#endif
		clipper.AddSubject(subjects);
		clipper.Execute(ClipType::Union, fillRule, unionSolution);

		// Then intersect with clip
		Paths64 solution;
		clipper.Clear();
		clipper.PreserveCollinear(preserveCollinear);
#ifdef USINGZ
		clipper.SetZCallback(UnionZCB);
#endif
		clipper.AddSubject(unionSolution);
		clipper.AddClip(clips);
		clipper.Execute(ClipType::Intersection, fillRule, solution);

		if (!solution.empty()) {
			outputPolylineCount = static_cast<int>(solution.size());

			*outputPolylines = (BoundaryPoint**)CoTaskMemAlloc(sizeof(BoundaryPoint*) * outputPolylineCount);
			*outputCounts = (int*)CoTaskMemAlloc(sizeof(int) * outputPolylineCount);

			for (int i = 0; i < outputPolylineCount; ++i) {
				int count = static_cast<int>(solution[i].size());
				(*outputCounts)[i] = count;

				(*outputPolylines)[i] = (BoundaryPoint*)CoTaskMemAlloc(sizeof(BoundaryPoint) * count);
				Path64ToBoundaryPoints(solution[i], (*outputPolylines)[i], count);
			}
		}
		else {
			outputPolylineCount = 0;
			*outputPolylines = nullptr;
			*outputCounts = nullptr;
		}
	}

} // namespace Clipper2Lib