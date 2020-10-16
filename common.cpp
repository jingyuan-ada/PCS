#include "common.h"


Common::Common()
{
}


Common::~Common()
{
}
//以K邻域最大距离计算每个点的密度,点数k/最大半径圆面积
float Common::computeLPD(const PointCloudT::ConstPtr cloud, int index, int k)
{
	float res = 0.0;
	vector<int> indices(k);
	vector<float> sqr_distances(k);
	search::KdTree<PointT> tree;
	tree.setInputCloud(cloud);
	if (tree.nearestKSearch(index, k, indices, sqr_distances) == k)
	{
		float dk = sqr_distances[k - 1];
		res = float(k) / static_cast<float> (M_PI)*dk*dk;
	}
	
	return res;
}
//k个邻近点距离和/k
float Common::computeMDK(const PointCloudT::ConstPtr cloud, int index, int k)
{
	float res = 0.0;
	vector<int> indices(k);
	vector<float> sqr_distances(k);
	search::KdTree<PointT> tree;
	tree.setInputCloud(cloud);
	if (tree.nearestKSearch(index, k, indices, sqr_distances) == k)
	{
		for (int i = 0; i < k; i++)
		{
			res += sqr_distances[i];
		}
	}

	return res/k;
}
//计算点云的平均距离，单点所占面积
float Common::computeCloudResolution(const PointCloudT::ConstPtr cloud, int k)
{
	float res = 0.0;
	int n_points = 0;
	int nres;
	vector<int> indices(k);
	vector<float> sqr_distances(k);
	search::KdTree<PointT> tree;
	tree.setInputCloud(cloud);
	for (size_t i = 0; i < cloud->size(); ++i)
	{
		if (!isfinite((*cloud)[i].x))
		{
			continue;
		}
		nres = tree.nearestKSearch(i, k, indices, sqr_distances);
		if (nres == k)
		{
			for (int i = 1; i < k; i++)
			{
				res += sqrt(sqr_distances[i]);
				++n_points;
			}
		}
	}
	if (n_points != 0)
	{
		res /= n_points;
	}
	return res;
}