#pragma once
#include <pcl/point_cloud.h>
#include <pcl/pcl_base.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>

using namespace std;
using namespace pcl;

typedef PointXYZ PointT;
typedef PointCloud<PointT> PointCloudT;

class Common
{
public:
	Common();
	~Common();
	float computeCloudResolution(const PointCloudT::ConstPtr cloud, int k);
	float computeLPD(const PointCloudT::ConstPtr cloud, int index, int k);//бшнд38
	float computeMDK(const PointCloudT::ConstPtr cloud, int index, int k);
private:
	

};

