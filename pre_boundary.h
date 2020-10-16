#pragma once
#ifndef PRE_BOUNDARY_H
#define PRE_BOUNDARY_H
#include <iostream>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/features/boundary.h>
#include <pcl/features/normal_3d.h>
#include <pcl/filters/normal_space.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <algorithm>

typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloudT;
typedef pcl::PointCloud<pcl::Normal> CloudNormal;
typedef pcl::PointCloud<pcl::Boundary> CloudBoundary;
class PreBoundary
{
public:
	PreBoundary() :
		normals_(new CloudNormal()),
		boundaries_()
		{ }
	void computeNormals(const PointCloudT::ConstPtr &cloud, float r);
	PointCloudT::Ptr computeCloudBoundary(const PointCloudT::ConstPtr &cloud);
	PointCloudT::Ptr computeCloudBoundary2(const PointCloudT::ConstPtr &cloud, int k, double m);
	PointCloudT::Ptr computeCloudBoundaryAGPN(const PointCloudT::ConstPtr &cloud, int k, double distance);
	bool isBoundary(const PointCloudT::ConstPtr &cloud, PointT q_point, std::vector<int> &indices, Eigen::Vector4f &u, const Eigen::Vector4f &v, float angle_threshold);
	~PreBoundary();

	CloudNormal::Ptr normals_;
	CloudBoundary boundaries_;

private:


};
#endif /*PRE_BOUNDARY_H*/
