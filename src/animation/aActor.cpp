#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

mat3 rotateAlign(vec3 u1, vec3 u2)
{
	vec3 axis = u1.Cross(u2);
	float dotProduct = u1 * u2;

	float angleRadians = acosf(dotProduct);

	const float sinA = sinf(angleRadians);
	const float cosA = cosf(angleRadians);
	const float invCosA = 1.0f - cosA;

	mat3 result( vec3((axis[0] * axis[0] * invCosA) + cosA,
		(axis[1] * axis[0] * invCosA) - (sinA * axis[2]),
		(axis[2] * axis[0] * invCosA) + (sinA * axis[1])),
		vec3((axis[0] * axis[1] * invCosA) + (sinA * axis[2]),
		(axis[1] * axis[1] * invCosA) + cosA,
		(axis[2] * axis[1] * invCosA) - (sinA * axis[0])),
		vec3((axis[0] * axis[2] * invCosA) - (sinA * axis[1]),
		(axis[1] * axis[2] * invCosA) + (sinA * axis[0]),
		(axis[2] * axis[2] * invCosA) + cosA)
	);

	return result;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	
	vec3 trans = m_pSkeleton->getRootNode()->getGlobalTranslation()  ;
	trans[2] = 0;
	m_Guide.setGlobalTranslation(trans);
	
	// 2.	Set the y component of the guide position to 0
	
	// 3.	Set the global rotation of the guide joint towards the guideTarget

	vec3 rd = m_Guide.getGlobalTranslation();

	m_Guide.setLocalRotation(m_Guide.getLocalRotation() * rotateAlign(rd, guideTargetPos));
	m_Guide.updateTransform();
	m_pSkeleton->update();
	
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);
	AJoint* root = m_pSkeleton->getRootNode();

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	m_pSkeleton->getRootNode()->setLocalTranslation(Max(leftHeight, rightHeight));
	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK 
	
	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal
		mat3 transform = root->getLocalRotation();
		transform[2] = leftNormal;
		root->setLocalRotation(root->getLocalRotation() * transform);
		root->updateTransform();
		
		
		;
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal
		;
	}
	m_pSkeleton->update();
}
