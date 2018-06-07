#include "kinematics_forward_solver.h"
#include <algorithm>
#include "math_utils.h"
#include "acclaim_skeleton.h"
#include "acclaim_motion.h"
#include "helper_forward_kinematics.h"
#include "kinematics_artic_idx.h"
#include "kinematics_pose.h"

using namespace math;

namespace kinematics {

// public func.

ForwardSolver::ForwardSolver()
    :skeleton_(nullptr),
    motion_(nullptr),
    artic_path_(new ArticIdxColl_t),
    helper_fk_(new helper::ForwardKinematics)
{
}

ForwardSolver::~ForwardSolver()
{
}

std::shared_ptr<acclaim::Skeleton> ForwardSolver::skeleton() const
{
    return skeleton_;
}

std::shared_ptr<acclaim::Motion> ForwardSolver::motion() const
{
    return motion_;
}

void ForwardSolver::set_skeleton(const std::shared_ptr<acclaim::Skeleton> &skeleton)
{
    skeleton_ = skeleton;
    helper_fk_->set_skeleton(skeleton_);
}

void ForwardSolver::set_motion(const std::shared_ptr<acclaim::Motion> &motion)
{
    motion_ = motion;
    //helper_fk_->set_skeleton(skeleton_);
	helper_fk_->set_motion(motion_);
}

void ForwardSolver::ConstructArticPath()
{
    helper_fk_->ConstructArticPath();
}

PoseColl_t ForwardSolver::ComputeSkeletonPose(const int32_t frame_idx)
{
    return this->ComputeSkeletonPose(motion_->joint_spatial_pos(frame_idx));
}

PoseColl_t ForwardSolver::ComputeSkeletonPose(const math::Vector6dColl_t &joint_spatial_pos)
{
	// TO DO

	PoseColl_t poseColl;

	Pose root;
	
	root.set_start_pos(joint_spatial_pos[0].linear_vector());
	root.set_end_pos(joint_spatial_pos[0].linear_vector());
	root.set_rotation(ComputeRotMatXyz(ToRadian(skeleton()->bone_ptr(0)->axis)));

	poseColl.push_back(root);

	for (int i = 1; i < skeleton()->bone_num(); i++)
	{
		acclaim::Bone currentBone = *skeleton()->bone_ptr(i);

		RotMat3d_t Rasf,Ramc,Ri0,Rroot;

		double tempIdentity[4][4];

		for (int x = 0; x < 4; x++)
		{
			for (int y = 0; y < 4; y++)
			{
				if (x == y)
				{
					tempIdentity[x][y] = 1;
				}
				else
				{
					tempIdentity[x][y] = 0;
				}
			}
		}

		Ri0 = ToRotMat(tempIdentity);

		Rroot = ComputeQuaternionXyz(
			ToRadian(joint_spatial_pos.at(0).angular_vector()).x(),
			ToRadian(joint_spatial_pos.at(0).angular_vector()).y(),
			ToRadian(joint_spatial_pos.at(0).angular_vector()).z());

		while (currentBone.idx != 0)
		{
			
			Rasf = ToRotMat(currentBone.rot_parent_current).transpose();

			Ramc = ComputeQuaternionXyz(
				ToRadian(joint_spatial_pos.at(currentBone.idx).angular_vector()).x(),
				ToRadian(joint_spatial_pos.at(currentBone.idx).angular_vector()).y(),
				ToRadian(joint_spatial_pos.at(currentBone.idx).angular_vector()).z());

			Ri0 = Rasf*Ramc*Ri0;

			currentBone = *currentBone.parent;
		}

		Ri0 = Rroot*Ri0;

		Vector3d_t Vi = skeleton()->bone_ptr(i)->dir*skeleton()->bone_ptr(i)->length;
		Vector3d_t Tim1;

		Tim1 = poseColl[skeleton()->bone_ptr(i)->parent->idx].end_pos();
		
		Vector3d_t Ti = Ri0*Vi + Tim1;

		Pose temp;
		temp.set_start_pos(Tim1);
		temp.set_end_pos(Ti);
		temp.set_rotation(ComputeRotMatXyz(ToRadian(skeleton()->bone_ptr(i)->axis)));

		poseColl.push_back(temp);
	}

	//return helper_fk_->ComputeSkeletonPose(joint_spatial_pos);
	return poseColl;
}

// protected func.

// private func.

} // namespace kinematics {
