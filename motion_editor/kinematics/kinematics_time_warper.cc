#include "kinematics_time_warper.h"
#include <utility>
#include "boost/numeric/conversion/cast.hpp"
#include "math_utils.h"

using namespace math;

namespace kinematics {

// public func.

TimeWarper::TimeWarper()
    :original_motion_sequence_(new math::SpatialTemporalVector6d_t),
    hard_constraint_coll_(new TimeWarpHardConstraintColl_t),
    time_step_(double{0.0}),
    min_time_step_(double{0.0}),
    max_time_step_(double{0.0})
{
}

TimeWarper::~TimeWarper()
{
}

double TimeWarper::time_step() const
{
    return time_step_;
}

double TimeWarper::min_time_step() const
{
    return min_time_step_;
}

double TimeWarper::max_time_step() const
{
    return max_time_step_;
}

void TimeWarper::Configure(
        const math::SpatialTemporalVector6d_t &original_motion_sequence,
        const double time_step,
        const double min_time_step,
        const double max_time_step
        )
{
    *original_motion_sequence_ = original_motion_sequence;
    time_step_ = time_step;
    min_time_step_ = min_time_step;
    max_time_step_ = max_time_step;
}

math::SpatialTemporalVector6d_t TimeWarper::ComputeWarpedMotion(
        const TimeWarpHardConstraintColl_t &hard_constraint_coll
        )
{
	// TO DO

	*hard_constraint_coll_ = hard_constraint_coll;
	SpatialTemporalVector6d_t answer = *original_motion_sequence_;

	int32_t originMidFrameidx = hard_constraint_coll_->at(2).frame_idx / hard_constraint_coll_->at(2).play_second * hard_constraint_coll_->at(1).play_second; // 150
	int32_t afterMidFrameidx = hard_constraint_coll_->at(1).frame_idx; // 160 or 140

	float index = (float)originMidFrameidx / (float)afterMidFrameidx;
	float index2 = (float)(hard_constraint_coll_->at(2).frame_idx - originMidFrameidx) / (float)(hard_constraint_coll_->at(2).frame_idx - afterMidFrameidx);

	for (int i = 0; i < original_motion_sequence_->temporal_size(); i++)
	{
		for (int j = 0; j < original_motion_sequence_->spatial_size(); j++)
		{
			if (i <= afterMidFrameidx)
			{
				Quaternion_t firstQ = ComputeQuaternionXyz(
					ToRadian(original_motion_sequence_->element(j, floorf(index*i)).angular_vector()).x(),
					ToRadian(original_motion_sequence_->element(j, floorf(index*i)).angular_vector()).y(),
					ToRadian(original_motion_sequence_->element(j, floorf(index*i)).angular_vector()).z());

				Quaternion_t secondQ = ComputeQuaternionXyz(
					ToRadian(original_motion_sequence_->element(j, ceilf(index*i)).angular_vector()).x(),
					ToRadian(original_motion_sequence_->element(j, ceilf(index*i)).angular_vector()).y(),
					ToRadian(original_motion_sequence_->element(j, ceilf(index*i)).angular_vector()).z());

				Quaternion_t::RealScalar slerpRatio = index*i - floorf(index*i);

				Vector6d_t temp = original_motion_sequence_->element(j, i);
				temp.set_angular_vector(ToDegree(ComputeEulerAngleXyz(ComputeRotMat(Slerp(firstQ, secondQ, slerpRatio)))));
				answer.set_element(j, i, temp);
			}
			else
			{
				Quaternion_t firstQ = ComputeQuaternionXyz(
					ToRadian(original_motion_sequence_->element(j, floorf(originMidFrameidx + index2*(i - afterMidFrameidx))).angular_vector()).x(),
					ToRadian(original_motion_sequence_->element(j, floorf(originMidFrameidx + index2*(i - afterMidFrameidx))).angular_vector()).y(),
					ToRadian(original_motion_sequence_->element(j, floorf(originMidFrameidx + index2*(i - afterMidFrameidx))).angular_vector()).z());
				Quaternion_t secondQ = ComputeQuaternionXyz(
					ToRadian(original_motion_sequence_->element(j, ceilf(originMidFrameidx + index2*(i - afterMidFrameidx))).angular_vector()).x(),
					ToRadian(original_motion_sequence_->element(j, ceilf(originMidFrameidx + index2*(i - afterMidFrameidx))).angular_vector()).y(),
					ToRadian(original_motion_sequence_->element(j, ceilf(originMidFrameidx + index2*(i - afterMidFrameidx))).angular_vector()).z());

				Quaternion_t::RealScalar slerpRatio = originMidFrameidx + index2*(i - afterMidFrameidx) - floorf(originMidFrameidx + index2*(i - afterMidFrameidx));

				Vector6d_t temp = original_motion_sequence_->element(j, i);
				temp.set_angular_vector(ToDegree(ComputeEulerAngleXyz(ComputeRotMat(Slerp(firstQ, secondQ, slerpRatio)))));
				answer.set_element(j, i, temp);
			}
		}
	}

	return answer;
}

// protected func.

// private func.

} // namespace kinematics {
