#include "sv_clustering.h"
#include "pre_boundary.h"
#include "boundary_VCCS.h"
#include <iostream>
#include <pcl/console/time.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/segmentation/lccp_segmentation.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/filter.h>
#include <pcl/visualization/pcl_visualizer.h>

using namespace std;
using namespace pcl;

void showCloud(const PointCloud<PointXYZL>::ConstPtr &cloud1, const PointCloud<PointXYZL>::ConstPtr &cloud2);
void showSingleCloud(const PointCloud<PointXYZ>::ConstPtr &cloud);

int main()
{
	PointCloud<PointXYZ>::Ptr cloud_input(new PointCloud<PointXYZ>());
	PreBoundary preBoundary;
	map<uint32_t, vector<BoundaryVCCS::BoundaryData>> cloud_boundary;
	//PointCloudT::Ptr cloud_boundary(new PointCloudT());
	float Rvoxel = 0.005f;//5
	float Rseed = 0.06f;//6
	//SVclustering svEst(Rvoxel,Rseed);
	io::loadPCDFile("pcdFile/test44.pcd", *cloud_input);
	//showSingleCloud(cloud_input);
	//对输入点云处理，滤除NAN点
	std::vector<int> mapping;
	removeNaNFromPointCloud(*cloud_input, *cloud_input, mapping);//转为dense
	console::TicToc tt;
	//double res = computeCloudResolution(cloud_input, 2);
	
	//*************get boundary points test****************

	//cloud_boundary = preBoundary.computeCloudBoundaryAGPN(cloud_input, 50, 0.01);
	//cloud_boundary = preBoundary.computeCloudBoundary2(cloud_input, 100, 7.2);
	preBoundary.computeNormals(cloud_input,Rvoxel);
	//cloud_boundary = preBoundary.computeCloudBoundary(cloud_input);
	//*************get supervoxel test**********************
	/*svEst.computeSupervoxel(cloud_input, preBoundary.normals_, preBoundary.boundaries_);
	std::vector<int> seed_indices;
	svEst.selectSeedVoxel(seed_indices);*/
	//svEst.createSVHelpers(seed_indices);
	BoundaryVCCS Bvccs(cloud_input,preBoundary.normals_);
	PointCloud<PointXYZL>::Ptr sv_labeled_cloud1(new PointCloud<PointXYZL>);
	PointCloud<PointXYZL>::Ptr sv_labeled_cloud2(new PointCloud<PointXYZL>);
	tt.tic();
	Bvccs.getVCCSs(Rvoxel, Rseed, sv_labeled_cloud1, sv_labeled_cloud2);
	//map<int, vector<Supervoxel<PointT>::Ptr>> regions1;
	map<int, vector<Supervoxel<PointT>::Ptr>> regions2;
	//Bvccs.svGrowing_c(regions1);
	
	Bvccs.svGrowing_rd(regions2);
	PointCloud<PointXYZL>::Ptr segment_cloud1(new PointCloud<PointXYZL>);
	segment_cloud1 = Bvccs.getPatchCloud();

	Bvccs.getRegionNeighbor();
	Bvccs.mergeSmallRegions();
	
	PointCloud<PointXYZL>::Ptr segment_cloud2(new PointCloud<PointXYZL>);
	segment_cloud2 = Bvccs.getPatchCloud();
	cout << "time:" << tt.toc() << "ms" << endl;
	//cout << cloud_boundary->size() << endl;
	
	showCloud(segment_cloud1, segment_cloud2);
	return 0;
}

//点云结果可视化
void showCloud(const PointCloud<PointXYZL>::ConstPtr &cloud1, const PointCloud<PointXYZL>::ConstPtr &cloud2)
{
	visualization::PCLVisualizer::Ptr viewer(new visualization::PCLVisualizer("cloud viewer"));
	int v1(0), v2(0);;
	viewer->createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer->setBackgroundColor(1.0, 1.0, 1.0, v1);
	viewer->addText("viewer1", 10, 10, "v1 text", v1);
	viewer->addPointCloud(cloud1, "cloud1", v1);
	viewer->createViewPort(0.5, 0.0, 1.0, 1.0, v2);
	viewer->setBackgroundColor(1.0, 1.0, 1.0, v2);
	viewer->addText("viewer2", 10, 10, "v2 text", v2);
	viewer->addPointCloud(cloud2, "cloud2", v2);
	viewer->setPointCloudRenderingProperties(visualization::PCL_VISUALIZER_POINT_SIZE, 2.0, "cloud1");
	viewer->setPointCloudRenderingProperties(visualization::PCL_VISUALIZER_POINT_SIZE, 2.0, "cloud2");
	while (!viewer->wasStopped())
	{
		//viewer->updatePointCloud(cloud, "cloud");
		//viewer->setPointCloudRenderingProperties(visualization::PCL_VISUALIZER_POINT_SIZE, 2.0, "cloud");
		viewer->spinOnce();
		//boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}
void showSingleCloud(const PointCloud<PointXYZ>::ConstPtr &cloud)
{
	visualization::PCLVisualizer::Ptr viewer(new visualization::PCLVisualizer("cloud viewer"));
	viewer->setBackgroundColor(1.0, 1.0, 1.0);
	visualization::PointCloudColorHandlerCustom<PointXYZ> single_color(cloud, 100, 100, 100);
	viewer->addPointCloud(cloud, single_color, "cloud1");
	//viewer->setPointCloudRenderingProperties(visualization::PCL_VISUALIZER_POINT_SIZE, 1.0, "cloud1");
	while (!viewer->wasStopped())
	{
		//viewer->updatePointCloud(cloud, "cloud");
		//viewer->setPointCloudRenderingProperties(visualization::PCL_VISUALIZER_POINT_SIZE, 2.0, "cloud");
		viewer->spinOnce();
		//boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}