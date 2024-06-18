
using UnityEngine;
using UnityEngine.Animations;
using System.Collections.Generic;
using System.Linq;

namespace QualisysRealTime.Unity{
public class AvatarMover : MonoBehaviour
{
    public Animator animator;

    private RTClient rtClient;
    private List<SixDOFBody> bodies;
    private HumanPoseHandler poseHandler;
    private HumanPose humanPose;

    private Dictionary<string, Queue<SixDOFBody>> bufferDict;
    private const int bufferSize = 5;

    void Start()
    {
        // Get the RTClient instance
        rtClient = RTClient.GetInstance();

        // Initialize the HumanPoseHandler
        if (animator.avatar != null && animator.avatar.isValid && animator.avatar.isHuman)
        {
            poseHandler = new HumanPoseHandler(animator.avatar, animator.transform);
            humanPose = new HumanPose();
        }
        else
        {
            Debug.LogError("Animator must have a valid humanoid avatar.");
        }

        // Set initial pose if needed
        SetInitialPose();

        // Initialize the buffer dictionary
        bufferDict = new Dictionary<string, Queue<SixDOFBody>>();
    }

    void Update()
    {
        if (poseHandler == null) return;

        // Fetch the updated SixDoF data from RTClient
        bodies = rtClient.Bodies;

        if (bodies != null && bodies.Count > 0)
        {
            // Update buffers with new data
            foreach (var body in bodies)
            {
                if (IsValidData(body))
                {
                    if (!bufferDict.ContainsKey(body.Name))
                    {
                        bufferDict[body.Name] = new Queue<SixDOFBody>();
                    }

                    var buffer = bufferDict[body.Name];
                    if (buffer.Count >= bufferSize)
                    {
                        buffer.Dequeue();
                    }
                    buffer.Enqueue(body);
                }
            }
        }

        // Update each body part with the corresponding buffered 6DoF data
        foreach (var bodyName in bufferDict.Keys)
        {
            HumanBodyBones bone = GetHumanBodyBone(bodyName);
            if (bone != HumanBodyBones.LastBone)
            {
                Transform bodyPart = animator.GetBoneTransform(bone);
                if (bodyPart != null)
                {
                    var bufferedBody = GetInterpolatedData(bodyName);
                    if (bufferedBody != null)
                    {
                        UpdateBodyPart(bodyPart, bufferedBody);
                    }
                }
            }
        }

        // Apply the updated pose to the avatar
        poseHandler.SetHumanPose(ref humanPose);
    }

    void UpdateBodyPart(Transform bodyPart, SixDOFBody data)
    {
        if (bodyPart != null && data != null)
        {
            bodyPart.localPosition = data.Position;
            bodyPart.localRotation = data.Rotation;
        }
    }

    void SetInitialPose()
    {
        if (poseHandler != null)
        {
            // Get the current pose
            poseHandler.GetHumanPose(ref humanPose);

            // Set the pose to T-pose or any default pose here if necessary
            // (For example, setting all muscle values to 0.0f for a T-pose)
            for (int i = 0; i < humanPose.muscles.Length; i++)
            {
                humanPose.muscles[i] = 0.0f;
            }

            // Apply the initial pose
            poseHandler.SetHumanPose(ref humanPose);
        }
    }

    HumanBodyBones GetHumanBodyBone(string bodyPartName)
    {
        switch (bodyPartName)
        {
            case "LeftUpperLeg":
                return HumanBodyBones.LeftUpperLeg;
            case "RightUpperLeg":
                return HumanBodyBones.RightUpperLeg;
            case "RightLowerLeg":
                return HumanBodyBones.RightLowerLeg;
            case "LeftLowerLeg":
                return HumanBodyBones.LeftLowerLeg;
            case "LeftHand":
                return HumanBodyBones.LeftHand;
            case "RightHand":
                return HumanBodyBones.RightHand;
            case "LeftFoot":
                return HumanBodyBones.LeftFoot;
            case "RightFoot":
                return HumanBodyBones.RightFoot;
            case "Head":
                return HumanBodyBones.Head;
            case "LeftLowerArm":
                return HumanBodyBones.LeftLowerArm;
            case "RightLowerArm":
                return HumanBodyBones.RightLowerArm;
            case "LeftUpperArm":
                return HumanBodyBones.LeftUpperArm;
            case "RightUpperArm":
                return HumanBodyBones.RightUpperArm;
            case "Hips":
                return HumanBodyBones.Hips; // Assuming torso refers to the hips
            case "Spine":
                return HumanBodyBones.Spine;
            default:
                return HumanBodyBones.LastBone; // Indicates an invalid bone
        }
    }

    SixDOFBody GetInterpolatedData(string bodyPartName)
    {
        if (bufferDict.ContainsKey(bodyPartName))
        {
            var buffer = bufferDict[bodyPartName];
            if (buffer.Count > 0)
            {
                if (buffer.Count == 1)
                {
                    return buffer.Peek(); // Only one data point available
                }
                else
                {
                    // Interpolate between the last two valid data points
                    SixDOFBody[] bufferArray = buffer.ToArray();
                    SixDOFBody last = bufferArray[bufferArray.Length - 1];
                    SixDOFBody secondLast = bufferArray[bufferArray.Length - 2];

                    if (IsValidData(secondLast) && IsValidData(last))
                    {
                        float t = 0.5f; // Midway interpolation
                        Vector3 interpolatedPosition = Vector3.Lerp(secondLast.Position, last.Position, t);
                        Quaternion interpolatedRotation = Quaternion.Slerp(secondLast.Rotation, last.Rotation, t);

                        return new SixDOFBody
                        {
                            Name = last.Name,
                            Position = interpolatedPosition,
                            Rotation = interpolatedRotation,
                            Color = last.Color
                        };
                    }
                    else
                    {
                        return IsValidData(last) ? last : secondLast;
                    }
                }
            }
        }
        return null;
    }

    bool IsValidData(SixDOFBody data)
    {
        return !float.IsNaN(data.Position.x) && !float.IsNaN(data.Position.y) && !float.IsNaN(data.Position.z)
            && !float.IsNaN(data.Rotation.x) && !float.IsNaN(data.Rotation.y) && !float.IsNaN(data.Rotation.z) && !float.IsNaN(data.Rotation.w);
    }
}
}
