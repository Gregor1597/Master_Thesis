using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using UnityEngine.UIElements;
using TMPro;
using UnityEditor;
using Unity.VRTemplate;
using QualisysRealTime.Unity;
using UnityEngine.InputSystem;
using Unity.VisualScripting;

public class Questionaire : MonoBehaviour
{

    private Dictionary<int,string> _questions = new Dictionary<int, string>(){
        //Ownership
        { 1, "It felt like the virtual body was my body." },
        { 2, "It felt like the virtual body parts were my body parts." },
        { 3, "The virtual body felt like a human body." },
        { 4, "It felt like the virtual body belonged to me." },
        //Agency
        { 5, "The movements of the virtual body felt like they were my movements." },
        { 6, "I felt like I was controlling the movements of the virtual body." },
        { 7, "I felt like I was causing the movements of the virtual body." },
        { 8, "The movements of the virtual body were in sync with my own movements." },
        //Change
        { 9, "I felt like the form or appearance of my own body had changed." },
        { 10, "I felt like the weight of my own body had changed." },
        { 11, "I felt like the size (height) of my own body had changed." },
        { 12, "I felt like the width of my own body had changed." },

    };

    public string ID = "test";
    private string _folder = "C:/Users/LauflaborVR2/Documents/Master_Gregor/data";

    public GameObject canvas;
    private Canvas myCanvas;

    public GameObject quest;

    private TextMeshProUGUI question;

    public GameObject val;
    private TextMeshProUGUI value;
    [HideInInspector]
    public int score; 

    [HideInInspector]   
    public bool nextCondition = false;
   //[HideInInspector]
    //public bool show = false; 
    private int clicked = 1; 

    public Experiment experiment;
    string condition = "Platzhalter";
    private int counter = 0;

    // Start is called before the first frame update
    void Start()
    {
        myCanvas =  canvas.GetComponent<Canvas>();
        question = quest.GetComponent<TextMeshProUGUI>();
        value = val.GetComponent<TextMeshProUGUI>();
        canvas.SetActive(false);
        condition = experiment.getCondition();
    }

    // Update is called once per frame
    void Update()
    {

                if(clicked == 1 && nextCondition== true ){
                        nextCondition = false;  
                         
                }
            
                 if(clicked > _questions.Count){
                    canvas.SetActive(false);
                    nextCondition = true; 
                    clicked = 1; 
                }
            

                if(Keyboard.current.escapeKey.isPressed){
                    //var rand = new System.Random();
                    //_questions = _questions.OrderBy(x => rand.Next()).ToDictionary(item => item.Key, item => item.Value);
                    canvas.SetActive(true);
                    //clicked = 1; 
                    question.text = _questions[clicked];

                }
                
                       
        
        
    }


public void OnButtonClicked (){
    //Debug.Log(_questions[clicked] + score);
    using(StreamWriter sw = new StreamWriter(_folder + "/" + ID + condition + "csv", append: true)){
        sw.WriteLine($"{clicked}, {score}");
    }
    clicked += 1;
    try{question.text = _questions[clicked];}
    catch(Exception ex){Debug.Log(ex);}
}

public void sliderScore(float sliderValue){
    score = ((int)sliderValue);
    value.text = score.ToString();
} 
    
}