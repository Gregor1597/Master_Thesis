using System;
using System.Collections;
using System.Collections.Generic;
using Unity.Tutorials.Core.Editor;
using Unity.VisualScripting;
using UnityEngine;
using QualisysRealTime.Unity; 
using Unity.XR.CoreUtils;
using UnityEditor;

public enum Identity{
        white,
        black,
        asian,
        MENA, 
        AIAN, 
        hispanic,
        NHPI
    };
public enum Gender {
        male, 
        female
    }; 
public class Experiment : MonoBehaviour
{
    
    List<string> conditions;
    [HideInInspector]
    public string condition;
    public Questionaire questionaire;
    private GameObject avatar;
    [HideInInspector]
    public GameObject tempGO; 
    private GameObject master; 
    public Identity identity; 
    public Gender gender;
    private List<GameObject> invisibles = new List<GameObject>(); 
    private int n; 
    public float smallScale = 0.8f; 
    public float largeScale = 1.2f; 
    public XROrigin origin; 
    private Transform CameraYOffset; 
    public Camera main; 
    private float startY;
    //public GameObject head;
    //public RTSkeleton stream; 
    // Start is called before the first frame update
    void Awake()
    {
        origin.MoveCameraToWorldLocation(new Vector3(0,0,0));
        createAvatar();
        conditions  = new List<string>(){
        
        "NoAvatar",
        "Normal",
        "Large",
        "Small"
        };
        foreach (Transform child in tempGO.transform){

            if(child.name != "master"||child.name != "Lights"){
                invisibles.Add(child.gameObject);
            }
            if (child.name == "master"){
                Debug.Log("Found master succesfully");
                master = child.gameObject;
            }
            if(child.name == "H_DDS_HighRes"){
               
                child.gameObject.AddComponent<BoxCollider>();
                child.gameObject.AddComponent<Rigidbody>();
                //child.gameObject.AddComponent<Collison_Handler>();
                //child.gameObject.GetComponent<BoxCollider>().isTrigger = true;
            }
        }
        CameraYOffset = origin.transform.GetChild(0); 
        startY = CameraYOffset.position.y; 
        var rand = new System.Random();
        n = rand.Next(conditions.Count);
        condition = conditions[n];
        //if a specific condition is needed
        //condition = "Normal"; 
        Debug.Log(condition);
        conditions.Remove(condition);
        changeAvatar(condition);
        
       

    }

    // Update is called once per frame
    void Update()
    {
        bool help = questionaire.nextCondition;
        /* muss noch implentiert werden
        if (stream.streaming == true && firstTime == true){
            origin.transform.position = head.transform.position - main.transform.position;
            firstTime = false; 
        }*/
        if(help == true){ 
            
            var rand = new System.Random();
            if(conditions.Count == 0){
                return; 
            }
            n = rand.Next(conditions.Count);
            condition = conditions[n];
            Debug.Log(condition);
            conditions.Remove(condition);
            changeAvatar(condition);

        }
       
        
    }

    /* nur eine idee, aber wahrscheinlich nicht
    private void scaleAvatar(float scale){
        Vector3 scaleTmp = object.transform.localScale;
        scaleTmp.x /= parent.localScale.x;
        scaleTmp.y /= parent.localScale.y;
        scaleTmp.z /= parent.localScale.z;
        object.transform.parent = parent;
        object.transform.localScale = scaleTmp;
    }*/

     private void changeAvatar(string condition){
        switch(condition){
            case "Normal":
                SetInviblesActive(true);
                master.transform.localScale = new Vector3(1, 1, 1);
                origin.transform.localScale = new Vector3(1,1,1); 
                break;
            case "Small":
                SetInviblesActive(true);
                master.transform.localScale = new Vector3(smallScale,smallScale,smallScale);
                origin.transform.localScale = new Vector3(smallScale,smallScale,smallScale); 
                
                break;
            case"Large":
                SetInviblesActive(true);
                master.transform.localScale = new Vector3(largeScale,largeScale,largeScale);
                origin.transform.localScale = new Vector3(largeScale,largeScale,largeScale); 
                
                break;
            case"NoAvatar":
                SetInviblesActive(false);
                origin.transform.localScale = new Vector3(1,1,1); 
                
                break;
            
        }
        origin.transform.localPosition = master.transform.position*master.transform.localScale.y;
    }
    public string getCondition(){
        return condition;
    }

    void SetInviblesActive(bool b){
        foreach(GameObject i in invisibles){
            i.SetActive(b);
        }
    }

    private void createAvatar(){
        Avatar destination = null;
        Debug.Log($"Idenity: {identity} , gender: {gender}");
        switch (identity){
            case Identity.white:
                if(gender == Gender.male)
                {
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/White_M_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/White_M_1_Casual.fbx",typeof(Avatar));
                }
                else{
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/White_F_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/White_F_1_Casual.fbx",typeof(Avatar));
                }
                break;
           case Identity.black:
                if(gender == Gender.male)
                {
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Black_M_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Black_M_1_Casual.fbx",typeof(Avatar));
                }
                else{
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Black_F_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Black_F_1_Casual.fbx",typeof(Avatar));
                }
                break;
            case Identity.asian:
                if(gender == Gender.male)
                {
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Asian_M_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Asian_M_1_Casual.fbx",typeof(Avatar));
                }
                else{
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Asian_F_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Asian_F_1_Casual.fbx",typeof(Avatar));
                }
                break;

            case Identity.MENA:
                if(gender == Gender.male)
                {
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/MENA_M_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/MENA_M_1_Casual.fbx",typeof(Avatar));
                }
                else{
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/MENA_F_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/MENA_F_1_Casual.fbx",typeof(Avatar));
                }
                break;
            case Identity.hispanic:
                if(gender == Gender.male)
                {
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Hispanic_M_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Hispanic_M_1_Casual.fbx",typeof(Avatar));
                }
                else{
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Hispanic_F_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/Hispanic_F_1_Casual.fbx",typeof(Avatar));
                }
                break;
            case Identity.NHPI:
                if(gender == Gender.male)
                {
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/NHPI_M_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/NHPI_M_1_Casual.fbx",typeof(Avatar));
                }
                else{
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/NHPI_F_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/NHPI_F_1_Casual.fbx",typeof(Avatar));
                }
                break;
            case Identity.AIAN:
                if(gender == Gender.male)
                {
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/AIAN_M_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/AIAN_M_1_Casual.fbx",typeof(Avatar));
                }
                else{
                    avatar = (GameObject)AssetDatabase.LoadAssetAtPath("Assets/Avatars/AIAN_F_1_Casual.fbx",typeof(GameObject));
                    destination = (Avatar)AssetDatabase.LoadAssetAtPath("Assets/Avatars/AIAN_F_1_Casual.fbx",typeof(Avatar));
                }
                break;
        }


        tempGO = Instantiate(avatar);
        ConvertToPrefabInstanceSettings convertToPrefabInstanceSettings = new ConvertToPrefabInstanceSettings();
        PrefabUtility.ConvertToPrefabInstance(tempGO, avatar, convertToPrefabInstanceSettings, InteractionMode.AutomatedAction);
        // add component to temp prefab instance
        tempGO.AddComponent<RTSkeleton>();
        // apply instance overrides to prefab
        //PrefabUtility.ApplyPrefabInstance(tempGO, InteractionMode.AutomatedAction);
        // destroy temp prefab instance
        //DestroyImmediate(tempGO);
        tempGO.GetComponent<RTSkeleton>().DestinationAvatar = destination;
        tempGO.GetComponent<RTSkeleton>().SkeletonName = "Q";//questionaire.ID;
        tempGO.transform.Rotate(0, 0, 0 );
        
    }
    
}

