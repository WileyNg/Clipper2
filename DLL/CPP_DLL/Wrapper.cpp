#pragma once
#define PINVOKE extern "C" __declspec(dllexport)

#include "clipper2/clipper.offset.h"

namespace Clipper2Lib {

	const double floating_point_tolerance = 1e-12;

	struct Point3d {
		double x;
		double y;
		double z;
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


	PINVOKE Point3d* OffsetPolylineOneSide(
		const Point3d* inputPoints,
		int inputCount,
		double delta,
		int& outputCount
	)
	{
		Path64 inputPath = ConvertToPath64(inputPoints, inputCount);

		ClipperOffset co;
		Paths64 subject;
		subject.push_back(inputPath);

		co.AddPaths(subject, JoinType::Miter, EndType::Butt);
		Paths64 solution;
		co.OffsetOpenPathOneSided(delta * 1e6, solution);
		
		if (solution.empty()) {
			outputCount = 0;
			return nullptr;
		}

		outputCount = static_cast<int>(solution[0].size());
		Point3d* outputPoints = new Point3d[outputCount];

		ConvertFromPath64(solution[0], outputPoints, outputCount);

		return outputPoints;
	}



	void ClipperOffset::OffsetOpenPathOneSided(double delta, Paths64& paths)
	{
		paths.clear();
		solution = &paths;
		solution->reserve(CalcSolutionCapacity());

		delta_ = delta;
		std::vector<Group>::iterator git;
		for (git = groups_.begin(); git != groups_.end(); ++git)
		{
			Group& group = *git;
			Paths64::const_iterator path_in_it = group.paths_in.cbegin();
			group_delta_ = std::abs(delta_);
			join_type_ = group.join_type;
			end_type_ = group.end_type;
			for (; path_in_it != group.paths_in.cend(); ++path_in_it) // DoGroupOffset(*git);
			{
				Path64::size_type pathLen = path_in_it->size();
				path_out.clear();
				const Path64& path = *path_in_it;
				BuildNormals(path);
				//OffsetOpenPath(group, path); 

				auto startPoint = path[0];
				auto startNormal = norms[0];
				path_out.push_back(Point64(startPoint.x + startNormal.x * delta, startPoint.y + startNormal.y * delta, startPoint.z));
				
				size_t highI = path.size() - 1;
				// offset the left side going forward
				for (Path64::size_type j = 1, k = 0; j < highI; k = j, ++j)
					OffsetPoint(group, path, j, k);


				auto endPoint = path[highI];
				auto endNormal = -norms[highI];
				path_out.push_back(Point64(endPoint.x + endNormal.x * delta, endPoint.y + endNormal.y * delta, endPoint.z));

				// reverse normals
				//for (size_t i = highI; i > 0; --i)
				//	norms[i] = PointD(-norms[i - 1].x, -norms[i - 1].y);
				//norms[0] = norms[highI];

				//Right side
				//for (size_t j = highI - 1, k = highI; j > 0; k = j, --j)
				//	OffsetPoint(group, path, j, k);
				solution->push_back(path_out);
			}
			
			if (!error_code_) 
				continue; // all OK
			solution->clear();
		}

		//if (!solution->size()) 
		//	return;

//		bool paths_reversed = CheckReverseOrientation();
//		//clean up self-intersections ...
//		Clipper64 c;
//		c.PreserveCollinear(false);
//		//the solution should retain the orientation of the input
//		c.ReverseSolution(reverse_solution_ != paths_reversed);
//#ifdef USINGZ
//		auto fp = std::bind(&ClipperOffset::ZCB, this, std::placeholders::_1,
//			std::placeholders::_2, std::placeholders::_3,
//			std::placeholders::_4, std::placeholders::_5);
//		c.SetZCallback(fp);
//#endif
//		c.AddSubject(*solution);
//		if (solution_tree)
//		{
//			if (paths_reversed)
//				c.Execute(ClipType::Union, FillRule::Negative, *solution_tree);
//			else
//				c.Execute(ClipType::Union, FillRule::Positive, *solution_tree);
//		}
//		else
//		{
//			if (paths_reversed)
//				c.Execute(ClipType::Union, FillRule::Negative, *solution);
//			else
//				c.Execute(ClipType::Union, FillRule::Positive, *solution);
//		}
	}
}