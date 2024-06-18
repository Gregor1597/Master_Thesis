using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using Unity.VisualScripting;
using QTMRealTimeSDK.Settings;
using System.Data.Common;
using System.Linq;
using System;
using Unity.XR.CoreUtils;

namespace QualisysRealTime.Unity
{
    public class RTselfmade : MonoBehaviour
    {
        private List<SixDOFBody> bodiesData;
        private RTClient rtClient;
        private GameObject markerRoot;
        private List<GameObject> bodies;

        public String SKname = "SR";

        private HumanPoseHandler mSourcePoseHandler;
        private HumanPoseHandler mDestiationPoseHandler;

        public Avatar DestinationAvatar;

        private Avatar mSourceAvatar;
        private HumanPose mHumanPose = new HumanPose();

        private Dictionary<string, string> mMecanimDictionary = new Dictionary<string, string>();

        private Dictionary<string, string> mMecanimDictionary2 = new Dictionary<string, string>();

        private Dictionary<string, string> BoneDictionary = new Dictionary<string, string>();

        private Dictionary<string,string> ParentDictionary = new Dictionary<string,string>();

        private Dictionary<uint, GameObject> QTMtoTransformDictionary;

        private GameObject mStreamedRootObject;

        private Skeleton ownSkeleton;//mQtmSkeletonCache;


        public bool visibleMarkers = true;

        [Range(0.001f, 1f)]
        public float markerScale = 0.05f;

        private bool streaming = false;

        

        void Start()
        {
            rtClient = RTClient.GetInstance();
            bodies = new List<GameObject>();
            markerRoot = gameObject;
            createMecanimDictionaries();
            createParentDictionary();
        }


        private void InitiateMarkers()
        {
            foreach (var marker in bodies)
            {
                Destroy(marker);
            }

            bodies.Clear();
            bodiesData = rtClient.Bodies;

            for (int i = 0; i < bodiesData.Count; i++)
            {
                GameObject newMarker = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                newMarker.name = bodiesData[i].Name;
                //Debug.Log("marker: " + newMarker.name + "position" + markerData[i].Position);
                newMarker.transform.parent = markerRoot.transform;
                newMarker.transform.localScale = Vector3.one * markerScale;
                newMarker.SetActive(false);
                bodies.Add(newMarker);
                
                
            }
            GameObject head = GameObject.CreatePrimitive(PrimitiveType.Sphere);
            head.name = "head";
            //Debug.Log("marker: " + head.name);
            head.transform.parent = markerRoot.transform;
            head.transform.localScale = Vector3.one * markerScale;
            head.SetActive(true);
            bodies.Add(head);
        }
        float timer = 0;
        bool time = true;
        // Update is called once per frame
        void Update()
        {
            if(time){
                timer += Time.deltaTime;
                    }
            if(time && timer > 5){ 
                time = false;
                if (rtClient == null) rtClient = RTClient.GetInstance();
                if (rtClient.GetStreamingStatus() && !streaming)
                {
                    InitiateMarkers();
                    streaming = true;
                }
                if (!rtClient.GetStreamingStatus() && streaming)
                {
                    streaming = false;
                    InitiateMarkers();
                }

                bodiesData = rtClient.Bodies;

                if (bodiesData == null || bodiesData.Count == 0)
                    return;

                if (bodies.Count != bodiesData.Count)
                {
                    InitiateMarkers();
                }

                
                for (int i = 0; i < bodiesData.Count; i++) 
                {
                    if (!float.IsNaN(bodiesData[i].Position.sqrMagnitude))
                    {
                        bodies[i].name = bodiesData[i].Name;
                        bodies[i].GetComponent<Renderer>().material.color = bodiesData[i].Color;
                        
                        if(bodies[i].name == "head") { bodies[i].transform.localPosition = bodiesData[5].Position + new Vector3(0,0.5f,0);}
                        else {bodies[i].transform.localPosition = bodiesData[i].Position;}
                        Debug.Log($" Marker :{bodies[i].name} with {bodiesData[i].Position} updated");
                        bodies[i].SetActive(true);
                        bodies[i].GetComponent<Renderer>().enabled = visibleMarkers;
                        bodies[i].transform.localScale = Vector3.one * markerScale;
                        
                    }
                    else
                    {
                        // Hide markers if we cant find them
                        bodies[i].SetActive(false);
                    }
                    
                    
                }
                if(ownSkeleton== null){
                    if(mStreamedRootObject != null)
                            GameObject.Destroy(mStreamedRootObject);
                        mStreamedRootObject = new GameObject(SKname);

                    BuildOwnSkeleton();
                    createMappingDictionary();
                    mStreamedRootObject.transform.SetParent(this.transform, false);
                    BuildMecanimAvatarFromQtmData();
                    return;
                }
                //Update all the game objects
                foreach (var segment in ownSkeleton.Segments)
                {
                    GameObject gameObject;
                    if (QTMtoTransformDictionary.TryGetValue(segment.Key, out gameObject))
                    {
                        gameObject.transform.localPosition = segment.Value.Position;
                        gameObject.transform.localRotation = segment.Value.Rotation;
                        Debug.Log($"GameObject {gameObject.name} updated at position {gameObject.transform.localPosition}");
                    }
                }

                if (mSourcePoseHandler != null && mDestiationPoseHandler != null)
                {
                    mSourcePoseHandler.GetHumanPose(ref mHumanPose);
                    mDestiationPoseHandler.SetHumanPose(ref mHumanPose);
                }
                }
        }
        /* private void BuildMecanimAvatarFromQtmData()

            {
                var humanBones =  new List<HumanBone>(ownSkeleton.Segments.Count);
                for (int index = 0; index < HumanTrait.BoneName.Length; index++)
                {
                    var humanBoneName = HumanTrait.BoneName[index];
                    
                    if (mMecanimDictionary.ContainsKey(humanBoneName))
                    {
                        
                        var bone = new HumanBone()
                        {
                            humanName = humanBoneName,
                            boneName = BoneDictionary[humanBoneName],
                        };
                        bone.limit.useDefaultValues = true;
                        humanBones.Add(bone);
                        //Debug.Log(bone.humanName + "   " + bone.boneName);
                    }
                }
                var skeletonBones = new List<SkeletonBone>(ownSkeleton.Segments.Count + 1);
                skeletonBones.Add(new SkeletonBone()
                {
                    name = this.SKname,
                    position = Vector3.zero,
                    // In QTM default poses are facing +Y which becomes -Z in Unity. 
                    // We rotate 180 degrees to match unity forward.
                    rotation = Quaternion.AngleAxis(180, Vector3.up), 
                    scale = Vector3.one,
                }); 

                foreach (var segment in ownSkeleton.Segments)
                {
                    skeletonBones.Add(new SkeletonBone()
                    {
                        name = segment.Value.Name,
                        position = segment.Value.TPosition,
                        rotation = segment.Value.TRotation,
                        scale = Vector3.one,
                    });
                }
            
                mSourceAvatar = AvatarBuilder.BuildHumanAvatar(mStreamedRootObject,
                    new HumanDescription()
                    {
                        human = humanBones.ToArray(),
                        skeleton = skeletonBones.ToArray(),
                        
                    }
                );
                if (mSourceAvatar.isValid == false || mSourceAvatar.isHuman == false)
                {
                    this.enabled = false;
                    return;
                }

                mSourcePoseHandler = new HumanPoseHandler(mSourceAvatar, mStreamedRootObject.transform);
                mDestiationPoseHandler = new HumanPoseHandler(DestinationAvatar, this.transform);
            }*/

        private void BuildMecanimAvatarFromQtmData()
        {
            var humanBones = new List<HumanBone>(HumanTrait.BoneName.Length);
            foreach (var humanBoneName in HumanTrait.BoneName)
            {
                if (BoneDictionary.TryGetValue(humanBoneName, out string boneName))
                {
                    if (!string.IsNullOrEmpty(boneName))
                    {
                        var bone = new HumanBone()
                        {
                            humanName = humanBoneName,
                            boneName = boneName,
                        };
                        bone.limit.useDefaultValues = true;
                        humanBones.Add(bone);
                        Debug.Log($"Added human bone: {bone.humanName} mapped to bone: {bone.boneName}");
                    }
                    else
                    {
                        Debug.LogError($"Bone name for human bone '{humanBoneName}' is empty or null in BoneDictionary.");
                    }
                }
                else
                {
                    Debug.LogWarning($"Human bone '{humanBoneName}' not found in BoneDictionary.");
                }
            }

            var skeletonBones = new List<SkeletonBone>(ownSkeleton.Segments.Count + 1);
            skeletonBones.Add(new SkeletonBone()
            {
                name = this.SKname,
                position = Vector3.zero,
                rotation = Quaternion.AngleAxis(180, Vector3.up),
                scale = Vector3.one,
            });
            Debug.Log($"Added Root bone: {this.SKname}");
            foreach (var segment in ownSkeleton.Segments)
            {
                skeletonBones.Add(new SkeletonBone()
                {
                    name = segment.Value.Name,
                    position = segment.Value.TPosition,
                    rotation = segment.Value.TRotation,
                    scale = Vector3.one,
                });
                Debug.Log($"Added skeleton bone: {segment.Value.Name} with position: {segment.Value.TPosition} and rotation: {segment.Value.TRotation}");
            }

            Debug.Log("Building Human Avatar...");
           // Log the hierarchy for debugging
            ValidateBoneHierarchy();
            mSourceAvatar = AvatarBuilder.BuildHumanAvatar(mStreamedRootObject,
                new HumanDescription()
                {
                    human = humanBones.ToArray(),
                    skeleton = skeletonBones.ToArray(),
                }
            );
            
            if (!mSourceAvatar.isValid || !mSourceAvatar.isHuman)
            {
                Debug.LogError("Built avatar is not valid or not human.");
                Debug.LogError($"IsValid: {mSourceAvatar.isValid}, IsHuman: {mSourceAvatar.isHuman}");
               
                
                this.enabled = false;
                return;
            }

            mSourcePoseHandler = new HumanPoseHandler(mSourceAvatar, mStreamedRootObject.transform);
            mDestiationPoseHandler = new HumanPoseHandler(DestinationAvatar, this.transform);

            
        }

        private void LogSkeletonHierarchy(GameObject root, string indent = "")
        {
            Debug.Log($"{indent}{root.name}");
            foreach (UnityEngine.Transform child in root.transform)
            {
                LogSkeletonHierarchy(child.gameObject, indent + "  ");
            }
        }

        private void ValidateBoneHierarchy()
        {
            // Check if all required bones are present in the hierarchy
            var requiredBones = new List<string>
            {
                "Hips", "Spine", "Chest", "UpperChest", "Neck", "Head",
                "LeftShoulder", "RightShoulder", "LeftUpperArm", "RightUpperArm",
                "LeftLowerArm", "RightLowerArm", "LeftHand", "RightHand",
                "LeftUpperLeg", "RightUpperLeg", "LeftLowerLeg", "RightLowerLeg",
                "LeftFoot", "RightFoot", "LeftToes", "RightToes"
            };

            foreach (var bone in requiredBones)
            {
                if (!BoneDictionary.ContainsKey(bone))
                {
                    Debug.LogError($"Missing required bone in BoneDictionary: {bone}");
                }
                else
                {
                    var transform = FindTransformRecursive(mStreamedRootObject.transform, BoneDictionary[bone]);
                    if (transform == null)
                    {
                        Debug.LogError($"Transform '{BoneDictionary[bone]}' for bone '{bone}' not found in hierarchy.");
                    }
                    else
                    {
                        Debug.Log($"Found transform '{BoneDictionary[bone]}' for bone '{bone}' in hierarchy.");
                    }
                }
            }
        }

        private UnityEngine.Transform FindTransformRecursive(UnityEngine.Transform parent, string name)
        {
            if (parent.name == name)
                return parent;

            foreach (UnityEngine.Transform child in parent)
            {
                var result = FindTransformRecursive(child, name);
                if (result != null)
                    return result;
            }

            return null;
        }


        private void createMecanimDictionaries(){
            mMecanimDictionary.Clear();
            mMecanimDictionary2.Clear();
            BoneDictionary.Clear();

            mMecanimDictionary.Add("Head", "head");
            //mMecanimDictionary.Add("Neck", "C7");

            //mMecanimDictionary.Add("Chest", "CLAV");
            //mMecanimDictionary.Add("UpperChest", "STRN");

            mMecanimDictionary.Add("Hips", "Hips"); // bei Testset aktuell RPSI und LPSI, aber lieber SAC
            mMecanimDictionary.Add("Spine", "Spine");

            mMecanimDictionary.Add("RightShoulder", "RSHO");
            mMecanimDictionary.Add("LeftShoulder", "LSHO");

            mMecanimDictionary.Add("RightUpperArm", "RightUpperArm"); //eigentlich ehr RUPA
            mMecanimDictionary.Add("LeftUpperArm", "LeftUpperArm");

            mMecanimDictionary.Add("RightLowerArm", "RightForeArm");
            mMecanimDictionary.Add("LeftLowerArm", "LeftLowerArm");

            mMecanimDictionary.Add("RightHand", "RightHand"); // oder RWRA
            mMecanimDictionary.Add("LeftHand", "LeftHand"); 

            mMecanimDictionary.Add("RightUpperLeg", "RightUpperLeg"); //alternativ RTHI
            mMecanimDictionary.Add("LeftUpperLeg", "LeftUpperLeg");
            
            mMecanimDictionary.Add("RightLowerLeg", "RightLowerLeg"); //alternativ RTIB
            mMecanimDictionary.Add("LeftLowerLeg", "LeftLowerLeg");

            mMecanimDictionary.Add("RightFoot", "RightFoot");
            mMecanimDictionary.Add("LeftFoot", "LeftFoot");

            //mMecanimDictionary.Add("RightToes", "RTOE");
           // mMecanimDictionary.Add("LeftToes", "LTOE");



            BoneDictionary.Add("Head", SKname + ":Head");
            BoneDictionary.Add("Neck", SKname + ":Neck");

            BoneDictionary.Add("Chest", SKname + ":Spine1");
            BoneDictionary.Add("UpperChest", SKname + ":Spine2");

            BoneDictionary.Add("Hips", SKname + ":Hips"); 
            BoneDictionary.Add("Spine", SKname + ":Spine");

            BoneDictionary.Add("RightShoulder", SKname + ":RightShoulder");
            BoneDictionary.Add("LeftShoulder", SKname + ":LeftShoulder");

            BoneDictionary.Add("RightUpperArm", SKname + ":RightArm"); //eigentlich ehr RUPA
            BoneDictionary.Add("LeftUpperArm", SKname + ":LeftArm");

            BoneDictionary.Add("RightLowerArm", SKname + ":RightForeArm");
            BoneDictionary.Add("LeftLowerArm", SKname + ":LeftForeArm");

            BoneDictionary.Add("RightHand", SKname + ":RightHand"); // oder RWRA
            BoneDictionary.Add("LeftHand", SKname + ":LeftHand"); 

            BoneDictionary.Add("RightUpperLeg", SKname + ":RightUpLeg"); //alternativ RTHI
            BoneDictionary.Add("LeftUpperLeg", SKname + ":LeftUpLeg");
            
            BoneDictionary.Add("RightLowerLeg", SKname + ":RightLeg"); //alternativ RTIB
            BoneDictionary.Add("LeftLowerLeg", SKname + "LeftLeg");

            BoneDictionary.Add("RightFoot", SKname + ":RightFoot");
            BoneDictionary.Add("LeftFoot", SKname + ":LeftFoot");

            //BoneDictionary.Add("RightToes", SKname + ":RightToeBase");
            //BoneDictionary.Add("LeftToes", SKname + ":LeftToeBase");

            mMecanimDictionary2 = mMecanimDictionary.ToDictionary(x => x.Value, x=> x.Key);
            //mMecanimDictionary.Add("Neck", "C7");

        }

        /*private void createMappingDictionary(){

            QTMtoTransformDictionary = new Dictionary<uint, GameObject>(ownSkeleton.Segments.Count);

            foreach(var segment in ownSkeleton.Segments){


                Debug.Log(segment.Value.Name + "");
                string n = mMecanimDictionary2[segment.Value.Name];

                string n2 = BoneDictionary[n];

                var gameObject = new GameObject(n2);
                
                gameObject.transform.parent = segment.Value.ParentId == 0 ? mStreamedRootObject.transform : QTMtoTransformDictionary[segment.Value.ParentId].transform;
                gameObject.transform.localPosition = segment.Value.TPosition;
                gameObject.transform.localRotation = segment.Value.TRotation;
                QTMtoTransformDictionary[segment.Value.Id] = gameObject;


            }
           
        }*/
        private void createMappingDictionary()
        {
            QTMtoTransformDictionary = new Dictionary<uint, GameObject>(ownSkeleton.Segments.Count);
            Debug.Log("Größe QTMDictionary: " + QTMtoTransformDictionary.Count);
            foreach (var segment in ownSkeleton.Segments){
                
                Debug.Log($"Segmente Key : Value --> {segment.Key} : {segment.Value.Name}");
            }
            
            foreach (var segment in ownSkeleton.Segments)
            {
                Debug.Log($"Segment: {segment.Value.Name}");
                string humanBoneName;
                if (mMecanimDictionary2.TryGetValue(segment.Value.Name, out humanBoneName))
                {
                    string boneName;
                    if (BoneDictionary.TryGetValue(humanBoneName, out boneName))
                    {
                        var gameObject = new GameObject(boneName);
                        gameObject.transform.parent = segment.Value.ParentId == 0 ? mStreamedRootObject.transform : QTMtoTransformDictionary[segment.Value.ParentId].transform;
                        gameObject.transform.localPosition = segment.Value.TPosition;
                        gameObject.transform.localRotation = segment.Value.TRotation;
                        QTMtoTransformDictionary[segment.Value.Id] = gameObject;
                        Debug.Log($"Created GameObject '{boneName}' for human bone '{humanBoneName}' with parent '{(segment.Value.ParentId == 0 ? mStreamedRootObject.name : QTMtoTransformDictionary[segment.Value.ParentId].name)} and position{gameObject.transform.localPosition}'");
                    }
                    else
                    {
                        Debug.LogError($"Bone name for '{humanBoneName}' not found in BoneDictionary.");
                    }
                }
                else
                {
                    Debug.LogError($"Segment name '{segment.Value.Name}' not found in mMecanimDictionary2.");
                }
            }
        }

        private void createParentDictionary(){
            ParentDictionary.Clear();

            //ParentDictionary.Add("Head","Neck");
            //ParentDictionary.Add("Neck","UpperChest");
           //ParentDictionary.Add("Head", "UpperChest");

            ParentDictionary.Add("Head","Spine");
            //ParentDictionary.Add("UpperChest","Chest");

            ParentDictionary.Add("Spine","Hips");
            ParentDictionary.Add("Hips","Reference");

            //ParentDictionary.Add("RightShoulder","UpperChest");
            //ParentDictionary.Add("LeftShoulder","UpperChest");

            ParentDictionary.Add("RightUpperArm","Spine");
            ParentDictionary.Add("LeftUpperArm","Spine");

            ParentDictionary.Add("RightLowerArm","RightUpperArm");
            ParentDictionary.Add("LeftLowerArm","LeftUpperArm");

            ParentDictionary.Add("RightHand","RightLowerArm");
            ParentDictionary.Add("LeftHand","LeftLowerArm");

            ParentDictionary.Add("LeftUpperLeg","Hips");
            ParentDictionary.Add("RightUpperLeg","Hips");

            ParentDictionary.Add("RightLowerLeg","RightUpperLeg");
            ParentDictionary.Add("LeftLowerLeg","LeftUpperLeg");

            ParentDictionary.Add("RightFoot","RightLowerLeg");
            ParentDictionary.Add("LeftFoot","LeftLowerLeg");

            //ParentDictionary.Add("RightToes","RightFoot");
            //ParentDictionary.Add("LeftToes","LeftFoot");

        }
    
        private void BuildOwnSkeleton (){

            if(ownSkeleton == null){

                ownSkeleton = new Skeleton();
            }
            else{
                ownSkeleton.Segments.Clear();
            }
            uint counter = 0;
            foreach (var m in bodies){
                if(mMecanimDictionary2.ContainsKey(m.name))
                {

                    Vector3 posi = m.transform.position;
                    Quaternion rot = m.transform.rotation;

                    Debug.Log($"Marker {m.name} with position {posi} and rotation {rot} added to dictionary");
                    string n = mMecanimDictionary2[m.name].ToString();
                    

                    Segment seg  = new Segment();
                    seg.Name = mMecanimDictionary[n];
                    seg.TPosition = posi;
                    seg.TRotation = rot;
                    seg.Id = counter;

                    //Debug.Log(n);
                

                    // eventuell hardcoden, dass dependencies richtig, vlt klärt sich dann das namefinding problem
                    ownSkeleton.Segments[counter] = seg;
                    counter++;
               
                }
            }
        
            foreach(var s in ownSkeleton.Segments){
                string parent = ParentDictionary[mMecanimDictionary2[s.Value.Name]];
                Debug.Log("parent:" + parent);
                if(parent == "Reference"){s.Value.ParentId = 0;}
                else{
                 for(uint i = 0; i < ownSkeleton.Segments.Count; i++)
                {
                    //Debug.Log("Check: " + mMecanimDictionary2[ownSkeleton.Segments[i].Name]);
                    if(mMecanimDictionary2[ownSkeleton.Segments[i].Name]==parent)
                    {
                        s.Value.ParentId = ownSkeleton.Segments[i].Id;
                        break;
                    }
                } 
                } 
            }

            
        }
    }
}