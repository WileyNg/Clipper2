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
		return Point64(pt.x + norm.x * delta, pt.y + norm.y * delta, pt.z);
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
	Paths64 SplitPathByConsecutiveTrueW(const Path64& path)
	{
		if (path.empty()) {
			return Paths64();
		}

		bool allTrue = true;

		Paths64 splitPaths;
		Path64 currentPath;
		Path64 leadingItems; // Items before first split marker
		bool foundFirstSplit = false;
		// Check if starts with single true and ends with single true
		bool startsWithSingleTrue = (path.size() > 1 && path[0].w && !path[1].w);
		bool endsWithSingleTrue = (path.size() > 1 && path[path.size() - 1].w && !path[path.size() - 2].w);
		bool specialCase = startsWithSingleTrue && endsWithSingleTrue;

		// Check if ends with double true (consecutive T,T at the end)
		bool endsWithDoubleTrue = (path.size() > 1 && path[path.size() - 1].w && path[path.size() - 2].w);

		for (size_t i = 0; i < path.size(); i++) {
			const Point64& pt = path[i];
			bool previousWasTrue = (i > 0) && path[i - 1].w;
			bool nextWasTrue = (i + 1 < path.size()) && path[i + 1].w;
			bool isIsolatedMiddleTrue;
			bool isConsecutiveTrue;
			// Check if current point starts a split marker sequence:
			// - Current is true AND (next is true OR it's isolated in middle)
			// - But NOT if previous was also true (we're in the middle of a sequence)
			bool isStartOfSplitMarker = false;
			if (pt.w && !previousWasTrue && i > 0) {
				// Current is true and previous was false
				// This is a split marker if:
				// 1. Next is also true (consecutive pair starts here)
				// 2. OR it's isolated in the middle (but NOT at the very end)
				isConsecutiveTrue = nextWasTrue;
				isIsolatedMiddleTrue = (i < path.size() - 1) && !nextWasTrue;

				isStartOfSplitMarker = isConsecutiveTrue || isIsolatedMiddleTrue;
			}
			// If we hit the start of a split marker
			if (isStartOfSplitMarker && !foundFirstSplit) {
				// First split found
				foundFirstSplit = true;
				if (!specialCase) {
					// Normal case: add split marker to current path, then save leading items
					currentPath.push_back(pt);
					leadingItems = currentPath;
					currentPath.clear();
				}
				else {
					// Special case: add the split marker to current path before saving
					currentPath.push_back(pt);
					if (!currentPath.empty()) {
						splitPaths.push_back(currentPath);
						currentPath.clear();
					}
				}
			}
			else if (isStartOfSplitMarker && foundFirstSplit) {
				// Save current path and start new one
				if (!currentPath.empty()) {
					if (!isIsolatedMiddleTrue)
						currentPath.push_back(pt);
					splitPaths.push_back(currentPath);
					currentPath.clear();
					if (isIsolatedMiddleTrue)
						currentPath.push_back(pt);
				}

				// If not the last split marker, skip the first T
				// If it is the last split marker, we already added it above
			}
			else {
				// Regular point (including the second T of consecutive pairs)
				currentPath.push_back(pt);
			}
		}
		// Handle the last path
		if (!currentPath.empty()) {
			if (!specialCase && !leadingItems.empty()) {
				// Check if we only have ONE split (not wrapping scenario)
				// This happens when there's only one split marker and the path doesn't wrap around
				bool hasIsolatedEndTrue = (path.size() > 0 && path[path.size() - 1].w &&
					(path.size() < 2 || !path[path.size() - 2].w));

				if (hasIsolatedEndTrue && splitPaths.empty()) {
					// Single split scenario - add paths separately without wrapping
					splitPaths.push_back(leadingItems);
					splitPaths.push_back(currentPath);
				}
				else {
					// Normal wrapping case: Prepend currentPath before leadingItems
					Path64 wrappedPath = currentPath;
					wrappedPath.insert(wrappedPath.end(), leadingItems.begin(), leadingItems.end());
					splitPaths.insert(splitPaths.begin(), wrappedPath);
				}
			}
			else {
				splitPaths.push_back(currentPath);
			}
		}
		return splitPaths;
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
	int getRequiredShift(Path64 path) {
		Point p1 = path[0];
		Point p2 = path[1];
		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;

		if (dx >= 0 && dy >= 0) return 0;  // Q1: x>0, y>0
		if (dx < 0 && dy >= 0)  return 1;  // Q2: x<0, y>0
		if (dx < 0 && dy < 0)   return 0;  // Q3: x<0, y<0
		return 1;                           // Q4: x>0, y<0
	}
	void ProcessOffsetSolution(
		Path64 path,
		const Paths64& solution,
		bool isClosedPath,
		bool smoothZ,
		bool simplifyPath,
		double epsilonForSimplifying,
		Point3d** outputPoints0,
		int& outputCount0,
		Point3d** outputPoints1,
		int& outputCount1
	)
	{
		// Initialize outputs
		*outputPoints0 = nullptr;
		*outputPoints1 = nullptr;
		outputCount0 = 0;
		outputCount1 = 0;

		Paths64 splitPaths;
		if (!isClosedPath&& solution.size() == 1 && solution[0].size() == 4) {

				int shift = getRequiredShift(path);
				// Split into two halves
				Path64 firstHalf;
				Path64 secondHalf;
				// First half: indices 1,2
				firstHalf.push_back(solution[0][(1 + shift) % 4]);
				firstHalf.push_back(solution[0][(2 + shift) % 4]);

				// Second half: indices 0,3 (with wrapping)
				secondHalf.push_back(solution[0][(0 + shift) % 4]);
				secondHalf.push_back(solution[0][(3 + shift) % 4]);

				splitPaths.push_back(firstHalf);
				splitPaths.push_back(secondHalf);
		}
		else {
			if (!isClosedPath)
				splitPaths = SplitPathByConsecutiveTrueW(solution[0]);
			else
			{
				Path64 firstPath = solution[0];
				firstPath.push_back(firstPath[0]);  // Append start point to close the path
				splitPaths.push_back(firstPath);

				if (solution.size() > 1) {
					Path64 lastPath = solution[solution.size() - 1];
					lastPath.push_back(lastPath[0]);  // Append start point to close the path
					splitPaths.push_back(lastPath);
				}
			}
		}
		if (splitPaths.empty()) {
			return;
		}

		// Process first path (splitPaths[0])
		if (splitPaths.size() == 1) {
			Path64 resultPath0 = splitPaths[0];
			if (smoothZ) {
				SmoothZValuesAcrossPlateaus(resultPath0);
			}
			if (simplifyPath) {
				resultPath0 = Clipper2Lib::SimplifyPath(resultPath0, epsilonForSimplifying, isClosedPath);
			}
			outputCount0 = static_cast<int>(resultPath0.size());
			*outputPoints0 = new Point3d[outputCount0];
			ConvertFromPath64(resultPath0, *outputPoints0, outputCount0);
		}
		// Process second path (splitPaths[1])
		else if (splitPaths.size() > 1)
		{
			Path64 resultPath0 = splitPaths[0];
			Path64 resultPath1 = (splitPaths.size() > 2) ? splitPaths[2] : splitPaths[1];

			// If there are 3 paths, append the last point from the second path to the first path
			if (splitPaths.size() > 2 && splitPaths[1].size() > 0) {
				resultPath0.push_back(splitPaths[1].back());
			}

			if (smoothZ) {
				SmoothZValuesAcrossPlateaus(resultPath0);
				SmoothZValuesAcrossPlateaus(resultPath1);
			}

			if (simplifyPath) {
				resultPath0 = Clipper2Lib::SimplifyPath(resultPath0, epsilonForSimplifying, isClosedPath);
				resultPath1 = Clipper2Lib::SimplifyPath(resultPath1, epsilonForSimplifying, isClosedPath);
			}

			outputCount0 = static_cast<int>(resultPath0.size());
			*outputPoints0 = new Point3d[outputCount0];

			outputCount1 = static_cast<int>(resultPath1.size());
			*outputPoints1 = new Point3d[outputCount1];

			ConvertFromPath64(resultPath0, *outputPoints0, outputCount0);
			ConvertFromPath64(resultPath1, *outputPoints1, outputCount1);
		}
	}

	// Exposed functions
	PINVOKE void OffsetPathConstant(
		const Point3d* inputPoints,
		int inputCount,
		double  delta,
		Point3d** outputPoints0,
		int& outputCount0,
		Point3d** outputPoints1,
		int& outputCount1,
		JoinType joinType = JoinType::Round,
		double epsilonForSimplifying = 10.0,
		bool simplifyBeforeOffset = true,
		bool smoothZ = true,
		bool simplifyPath = true
	)
	{
		// Initialize outputs
		*outputPoints0 = nullptr;
		*outputPoints1 = nullptr;
		outputCount0 = 0;
		outputCount1 = 0;

		bool isClosedPath = IsPathClosed(inputPoints, inputCount);
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);

		if (simplifyBeforeOffset) {
			inputPath = Clipper2Lib::SimplifyPath(inputPath, epsilonForSimplifying, isClosedPath);
		}

		ClipperOffset co;
		co.IsClosed = isClosedPath;
		Paths64 subject;
		subject.push_back(inputPath);
 		co.AddPaths(subject, joinType, isClosedPath?  EndType::Polygon: EndType::Butt);

		Paths64 solution;

		//co.SetDeltaCallback([delta](const Path64& path,
		//	const PathD& path_norms, size_t curr_idx, size_t prev_idx)
		//	{
		//		// gradually scale down the offset to a minimum of 25% of delta
		//		double high = static_cast<double>(path.size() - 1) * 1.25;
		//		return (curr_idx-high) / high * delta*1e6;
		//	});

		co.OffsetPathConstant(delta * 1e6, solution);

		if (isClosedPath)
		{
			Paths64 solution2;
			co.OffsetPathConstant(delta * 1e6*-1, solution2);
			if (solution2.size() > 0)
				solution.push_back(solution2[0]);
		}

		if (solution.empty()) {
			return;
		}

		ProcessOffsetSolution(inputPath, solution, isClosedPath, smoothZ, simplifyPath,
			epsilonForSimplifying, outputPoints0, outputCount0,
			outputPoints1, outputCount1);
	}

	// Exposed functions
	PINVOKE void OffsetPathVariable(
		const Point3d* inputPoints,
		int inputCount,
		double* deltas,
		Point3d** outputPoints0,
		int& outputCount0,
		Point3d** outputPoints1,
		int& outputCount1,
		JoinType joinType = JoinType::Round,
		double epsilonForSimplifying = 10.0,
		bool simplifyBeforeOffset = true,
		bool smoothZ = true,
		bool simplifyPath = true
	)
	{
		// Initialize outputs
		*outputPoints0 = nullptr;
		*outputPoints1 = nullptr;
		outputCount0 = 0;
		outputCount1 = 0;

		bool isClosedPath = IsPathClosed(inputPoints, inputCount);
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);

		if (simplifyBeforeOffset) {
			inputPath = Clipper2Lib::SimplifyPath(inputPath, epsilonForSimplifying, isClosedPath);
		}

		ClipperOffset co;
		Paths64 subject;
		Paths64 splitPaths;
		subject.push_back(inputPath);
		co.AddPaths(subject, joinType, isClosedPath ?  EndType::Polygon: EndType::Round);

		Paths64 solution;

		co.SetDeltaCallback([deltas](const Path64& path,
			const PathD& path_norms, size_t curr_idx, size_t prev_idx)
			{
				return  deltas[curr_idx] * 1e6;
			});

		co.OffsetPathVariable(solution);
		if (isClosedPath)
		{
			Paths64 solution2;
			co.SetDeltaCallback([deltas](const Path64& path,
				const PathD& path_norms, size_t curr_idx, size_t prev_idx)
				{
					return  deltas[curr_idx] * 1e6*-1;
				});
			co.OffsetPathVariable(solution2);
			solution.push_back(solution2[0]);
		}

		if (solution.empty()) {
			return;
		}
		ProcessOffsetSolution(inputPath, solution, isClosedPath, smoothZ, simplifyPath,
			epsilonForSimplifying, outputPoints0, outputCount0,
			outputPoints1, outputCount1);
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
			Paths64{inputPath}, // expects Paths64
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
		} else {
			outputCount0 = 0;
			*outputPoints0 = nullptr;
		}
	}
	PINVOKE void InflateVariable(
		const Point3d* inputPoints,
		int inputCount,
		double*  deltas,
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
				return  deltas[curr_idx] * 1e6  ;
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
		} else {
			outputCount0 = 0;
			*outputPoints0 = nullptr;
		}
	}
	extern "C" __declspec(dllexport)
		void FreeMemory(Point3d* ptr)
	{
		delete[] ptr;
	}


	void ClipperOffset::OffsetPathConstant(double delta, Paths64& paths)
	{
		//Begins copy of ExecuteInternal(double delta) ; 

		paths.clear();
		solution = &paths;
		solution->reserve(CalcSolutionCapacity());
#
		if (std::abs(delta) < 0.5) // ie: offset is insignificant
		{
			Paths64::size_type sol_size = 0;
			for (const Group& group : groups_) sol_size += group.paths_in.size();
			solution->reserve(sol_size);
			for (const Group& group : groups_)
				copy(group.paths_in.begin(), group.paths_in.end(), back_inserter(*solution));
		}
		else
		{
			temp_lim_ = (miter_limit_ <= 1) ?
				2.0 :
				2.0 / (miter_limit_ * miter_limit_);

			delta_ = delta;
			std::vector<Group>::iterator git;
			for (git = groups_.begin(); git != groups_.end(); ++git)
			{
				Group& group = *git;
				Paths64::const_iterator path_in_it = group.paths_in.cbegin();

				// Begins copy of DoGroupOffset
				if (group.end_type == EndType::Polygon)
				{
					// a straight path (2 points) can now also be 'polygon' offset
					// where the ends will be treated as (180 deg.) joins
					if (!group.lowest_path_idx.has_value()) delta_ = std::abs(delta_);
					group_delta_ = (group.is_reversed) ? -delta_ : delta_;
				}
				else
					group_delta_ = std::abs(delta_);// *0.5;

				double abs_delta = std::fabs(group_delta_);


				join_type_ = group.join_type;
				end_type_ = group.end_type;
				if (group.join_type == JoinType::Round || group.end_type == EndType::Round)
				{
					double arcTol = (arc_tolerance_ > floating_point_tolerance ?
						std::min(abs_delta, arc_tolerance_) :
						std::log10(2 + abs_delta) * default_arc_tolerance);

					double steps_per_360 = std::min(PI / std::acos(1 - arcTol / abs_delta), abs_delta * PI);
					step_sin_ = std::sin(2 * PI / steps_per_360);
					step_cos_ = std::cos(2 * PI / steps_per_360);
					if (group_delta_ < 0.0) step_sin_ = -step_sin_;
					steps_per_rad_ = steps_per_360 / (2 * PI);
				}

				for (; path_in_it != group.paths_in.cend(); ++path_in_it)
				{
					Path64::size_type pathLen = path_in_it->size();

					if ((pathLen == 2) && (group.end_type == EndType::Joined))
						end_type_ = (group.join_type == JoinType::Round) ?
						EndType::Round :
						EndType::Square;

					const Path64& path = *path_in_it;
					BuildNormals(path);

					path_out.clear();

					// Extracted method call
					if(IsClosed)
						OffsetPolygon(group, path);
					else
						OffsetOpenPathWithZ(group, path);
					solution->push_back(path_out);
				}
				//Begins copy of rest of ExecuteInternal(double delta) ; 

				if (!solution->size()) return;

				bool paths_reversed = CheckReverseOrientation();
				//clean up self-intersections ...
				Clipper64 c;
				c.PreserveCollinear(false);
				//the solution should retain the orientation of the input
				c.ReverseSolution(reverse_solution_ != paths_reversed);
#ifdef USINGZ
				auto fp = std::bind(&ClipperOffset::ZCB, this, std::placeholders::_1,
					std::placeholders::_2, std::placeholders::_3,
					std::placeholders::_4, std::placeholders::_5);
				c.SetZCallback(fp);
#endif
				c.AddSubject(*solution);
				if (solution_tree)
				{
					if (paths_reversed)
						c.Execute(ClipType::Union, FillRule::Negative, *solution_tree);
					else
						c.Execute(ClipType::Union, FillRule::Positive, *solution_tree);
				}
				else
				{
					if (paths_reversed)
						c.Execute(ClipType::Union, FillRule::Negative, *solution);
					else
						c.Execute(ClipType::Union, FillRule::Positive, *solution);
				}
			}

		}
	}

	void ClipperOffset::OffsetOpenPathWithZ(Group& group, const Path64& path)
	{
		auto startPoint = path[0];
		auto startNormal = norms[0];
		double abs_delta = std::abs(group_delta_);
		auto pt1 = Point64(startPoint.x - abs_delta * startNormal.x, startPoint.y - abs_delta * startNormal.y, startPoint.z, true);
		auto pt2 = Point64(startPoint.x + abs_delta * startNormal.x, startPoint.y + abs_delta * startNormal.y, startPoint.z, true);

		path_out.push_back(Point64(pt1));
		path_out.push_back(Point64(pt2));

		size_t highI = path.size() - 1;

		// offset the left side going forward
		// k is current, j is next
		for (Path64::size_type j = 1, k = 0; j < highI; k = j, ++j)
			OffsetPointWithZ(group, path, j, k);

		// reverse normals
		for (size_t i = highI; i > 0; --i)
			norms[i] = PointD(-norms[i - 1].x, -norms[i - 1].y);
		norms[0] = norms[highI];

		// do the line end cap
		if (deltaCallback64_)
			group_delta_ = deltaCallback64_(path, norms, highI, highI);

		if (std::fabs(group_delta_) <= floating_point_tolerance)
			path_out.push_back(path[highI]);
		else
		{
			double abs_delta = std::abs(group_delta_);

			auto pt1 = Point64(path[highI].x - abs_delta * norms[highI].x, path[highI].y - abs_delta * norms[highI].y, path[highI].z, true);
			auto pt2 = Point64(path[highI].x + abs_delta * norms[highI].x, path[highI].y + abs_delta * norms[highI].y, path[highI].z, true);
			path_out.push_back(Point64(pt1));
			path_out.push_back(Point64(pt2));
		}

		// offset the right side going backward
		for (size_t j = highI - 1, k = highI; j > 0; k = j, --j)
			OffsetPointWithZ(group, path, j, k);
	}

	void ClipperOffset::OffsetPathVariable(Paths64& paths)
	{
		//Begins copy of ExecuteInternal(double delta) ; 

		paths.clear();
		solution = &paths;
		solution->reserve(CalcSolutionCapacity());


		temp_lim_ = (miter_limit_ <= 1) ?
			2.0 :
			2.0 / (miter_limit_ * miter_limit_);


		std::vector<Group>::iterator git;
		for (git = groups_.begin(); git != groups_.end(); ++git)
		{
			// Begins copy of DoGroupOffset(*git);
			Group& group = *git;
			Paths64::const_iterator path_in_it = group.paths_in.cbegin();
			group_delta_ = std::abs(delta_);
			join_type_ = group.join_type;
			end_type_ = group.end_type;

			for (; path_in_it != group.paths_in.cend(); ++path_in_it)
			{
				Path64::size_type pathLen = path_in_it->size();
				path_out.clear();
				if ((pathLen == 2) && (group.end_type == EndType::Joined))
					end_type_ = (group.join_type == JoinType::Round) ?
					EndType::Round :
					EndType::Square;

				const Path64& path = *path_in_it;
				BuildNormals(path);
				if (deltaCallback64_)
					group_delta_ = deltaCallback64_(path, norms, 0, 0);


				//Begins copy of OffsetOpenPath(group, *path_in_it);

				if (IsClosed)
					OffsetPolygon(group, path);
				else
					OffsetOpenPathWithZ(group, path);
				solution->push_back(path_out);
			}
			//Begins copy of rest of ExecuteInternal(double delta) ; 

			if (!solution->size()) return;

			bool paths_reversed = CheckReverseOrientation();
			//clean up self-intersections ...
			Clipper64 c;
			c.PreserveCollinear(false);
			//the solution should retain the orientation of the input
			c.ReverseSolution(reverse_solution_ != paths_reversed);

			auto fp = std::bind(&ClipperOffset::ZCB, this, std::placeholders::_1,
				std::placeholders::_2, std::placeholders::_3,
				std::placeholders::_4, std::placeholders::_5);
			c.SetZCallback(fp);

			c.AddSubject(*solution);
			if (solution_tree)
			{
				if (paths_reversed)
					c.Execute(ClipType::Union, FillRule::Negative, *solution_tree);
				else
					c.Execute(ClipType::Union, FillRule::Positive, *solution_tree);
			}
			else
			{
				if (paths_reversed)
					c.Execute(ClipType::Union, FillRule::Negative, *solution);
				else
					c.Execute(ClipType::Union, FillRule::Positive, *solution);
			}
		}

	}
	void ClipperOffset::OffsetPointWithZ(Group& group, const Path64& path, size_t j, size_t k)
	{
		//j is current, k is previous

		// Let A = change in angle where edges join
		// A == 0: ie no change in angle (flat join)
		// A == PI: edges 'spike'
		// sin(A) < 0: right turning
		// cos(A) < 0: change in angle is more than 90 degree

		if (path[j] == path[k]) return;

		double sin_a = CrossProduct(norms[j], norms[k]);
		double cos_a = DotProduct(norms[j], norms[k]);
		if (sin_a > 1.0) sin_a = 1.0;
		else if (sin_a < -1.0) sin_a = -1.0;

		if (deltaCallback64_) {
			group_delta_ = deltaCallback64_(path, norms, j, k);
			if (group.is_reversed)
				group_delta_ = -group_delta_;
		}
		if (std::fabs(group_delta_) <= floating_point_tolerance)
		{
			path_out.push_back(path[j]);
			return;
		}

		if (cos_a > -0.999 && (sin_a * group_delta_ < 0)) // test for concavity first (#593)
		{

			// is concave
			// by far the simplest way to construct concave joins, especially those joining very 
			// short segments, is to insert 3 points that produce negative regions. These regions 
			// will be removed later by the finishing union operation. This is also the best way 
			// to ensure that path reversals (ie over-shrunk paths) are removed.
			Point64 tmp = GetPerpendic(path[j], norms[k], group_delta_);
			tmp.z = path[j].z;
			path_out.push_back(tmp);
			// when the angle is almost flat (cos_a ~= 1), it's safe to skip this middle point
			if (cos_a < 0.999) path_out.push_back(path[j]); // (#405, #873)

			Point64 tmp2 = GetPerpendic(path[j], norms[j], group_delta_);
			tmp2.z = path[j].z;
			path_out.push_back(tmp2);
		}
		else if (cos_a > 0.999 && join_type_ != JoinType::Round)
		{
			// almost straight - less than 2.5 degree (#424, #482, #526 & #724)
			DoMiter(path, j, k, cos_a);
		}
		else if (join_type_ == JoinType::Miter)
		{
			// miter unless the angle is sufficiently acute to exceed ML
			if (cos_a > temp_lim_ - 1) DoMiter(path, j, k, cos_a);
			else DoSquare(path, j, k);
		}
		else if (join_type_ == JoinType::Round)
			DoRound(path, j, k, std::atan2(sin_a, cos_a));
		else if (join_type_ == JoinType::Bevel)
			DoBevel(path, j, k);
		else
			DoSquare(path, j, k);
	}

#include <vector>
#include <cmath>


	

}