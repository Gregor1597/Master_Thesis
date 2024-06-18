using System;
using System.Collections;
using System.Collections.Generic;
using Unity.Tutorials.Core.Editor;
using Unity.VisualScripting;
using UnityEngine;


public class Experiment : MonoBehaviour
{

    List<string> conditions;
    [HideInInspector]
    public string condition;

    public Questionaire questionaire;
    public GameObject avatar; 
    private int n; 
    public float smallScale = 0.8f; 
    public float largeScale = 1.2f; 
    // Start is called before the first frame update
    void Start()
    {
        conditions  = new List<string>(){
        "NoAvatar",
        "Normal",
        "Large",
        "Small"
        };
        var rand = new System.Random();
        n = rand.Next(conditions.Count);
        condition = conditions[n];
        conditions.Remove(condition);
        changeAvatar(condition);

    }

    // Update is called once per frame
    void Update()
    {
        bool help = questionaire.nextCondition;
        
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
                avatar.SetActive(true);
                avatar.transform.localScale = new Vector3(1, 1, 1);
                 break;
            case "Small":
                avatar.SetActive(true);
                avatar.transform.localScale = new Vector3(smallScale,smallScale,smallScale);
                break;
            case"Large":
                avatar.SetActive(true);
                avatar.transform.localScale = new Vector3(largeScale,largeScale,largeScale);
                break;
            case"NoAvatar": 
                avatar.SetActive(false);
                break;
            
        }
    }
    public string getCondition(){
        return condition;
    }
}
