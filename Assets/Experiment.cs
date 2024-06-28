using System;
using System.Collections;
using System.Collections.Generic;
using Unity.Tutorials.Core.Editor;
using Unity.VisualScripting;
using UnityEngine;
using QualisysRealTime.Unity; 
using Unity.XR.CoreUtils; 


public class Experiment : MonoBehaviour
{
    
    List<string> conditions;
    [HideInInspector]
    public string condition;
    public bool male; 
    public Questionaire questionaire;
    public GameObject Mavatar; 
    public GameObject FAvatar; 
    private GameObject avatar;
    public GameObject master; 
    private List<GameObject> invisibles = new List<GameObject>(); 
    private int n; 
    public float smallScale = 0.8f; 
    public float largeScale = 1.2f; 
    public XROrigin origin; 
    private Transform CameraYOffset; 
    public Camera main; 
    private float startY;
    public GameObject head;
    public RTSkeleton stream; 
    private bool firstTime;
    // Start is called before the first frame update
    void Start()
    {
        origin.MoveCameraToWorldLocation(new Vector3(0,0,0));
        avatar = (male) ? Mavatar : FAvatar;
        conditions  = new List<string>(){
        "NoAvatar",
        "Normal",
        "Large",
        "Small"
        };
        foreach (Transform child in avatar.transform){
            if(child.tag =="Face"){
                invisibles.Add(child.gameObject);
            }
        }
        CameraYOffset = origin.transform.GetChild(0); 
        startY = CameraYOffset.position.y; 
        var rand = new System.Random();
        n = rand.Next(conditions.Count);
        condition = conditions[n];
        Debug.Log(condition);
        conditions.Remove(condition);
        changeAvatar(condition);
        Debug.Log(head.transform.position);
        firstTime = false;
       

    }

    // Update is called once per frame
    void Update()
    {
        bool help = questionaire.nextCondition;
         if (stream.streaming == true && firstTime == true){
            origin.transform.position = head.transform.position - main.transform.position;
            firstTime = false; 
        }
        if(help == true){
            
            var rand = new System.Random();
            if(conditions.Count == 0){
                return; 
            }
            n = rand.Next(conditions.Count);
            condition = conditions[n];
            conditions.Remove(condition);
            changeAvatar(condition);

        }
       
        
    }

    private void changeAvatar(string condition){
        switch(condition){
            case "Normal":
                SetInviblesActive(true);
                master.transform.localScale = new Vector3(1, 1, 1);
                //CameraYOffset.position = new Vector3(CameraYOffset.position.x, startY, CameraYOffset.position.z);
                break;
            case "Small":
                SetInviblesActive(true);
                master.transform.localScale = new Vector3(smallScale,smallScale,smallScale);
                //CameraYOffset.position = new Vector3(CameraYOffset.position.x, startY + smallScale, CameraYOffset.position.z);
                break;
            case"Large":
                SetInviblesActive(true);
                master.transform.localScale = new Vector3(largeScale,largeScale,largeScale);
                //CameraYOffset.position = new Vector3(CameraYOffset.position.x, startY * largeScale, CameraYOffset.position.z);
                break;
            case"NoAvatar":
                SetInviblesActive(false);
                //CameraYOffset.position = new Vector3(CameraYOffset.position.x, startY, CameraYOffset.position.z);
                break;
            
        }
        /*if (stream.streaming == true){
            origin.transform.position = head.transform.position - main.transform.localPosition;
        }*/ //funktioniert nicht 
    }
    public string getCondition(){
        return condition;
    }

    void SetInviblesActive(bool b){
        foreach(GameObject i in invisibles){
            i.SetActive(b);
        }
    }
}
