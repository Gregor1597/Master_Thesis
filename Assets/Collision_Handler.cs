using System.Collections;
using System.Collections.Generic;
using QTMRealTimeSDK.Settings;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements; 
using TMPro;
using UnityEngine.InputSystem;

public class Collision_Handler : MonoBehaviour
{

    [HideInInspector]
    public int counter; 
    
    public GameObject text_mesh; 
    public GameObject black_mesh;
    private TextMeshProUGUI  text;

    private float waitingTime = 1f;
    private bool black_screen_active = false;
    // Start is called before the first frame update
    void Start()
    {
        counter = 0;
        black_mesh.SetActive(black_screen_active);
        text = text_mesh.GetComponent<TextMeshProUGUI>();
        text_mesh.SetActive(false);
    }

    // Update is called once per frame
    void Update()
    {
        if(Keyboard.current.spaceKey.isPressed ){
            black_screen_active = !black_screen_active;
            black_mesh.SetActive(black_screen_active);
        }
    }

    private void OnTriggerEnter(Collider collision){

        counter += 1; 
        StartCoroutine(ShowAndWait(waitingTime));
        StartCoroutine(DisableCollider(5*waitingTime));
        Debug.Log("Triggered,, counter = " + counter );
        text.text = counter.ToString();
    }
    public int getCounter(){
        return counter;
    }
    public void setCounter(int c){
        counter = c;
    }
    private IEnumerator ShowAndWait(float waitTime){
        text_mesh.SetActive(true);
        yield return new WaitForSeconds(waitTime);
        text_mesh.SetActive(false);

    }
    private IEnumerator DisableCollider(float waitTime){
           this.gameObject.SetActive(false);
           yield return new WaitForSeconds(waitTime);
            this.gameObject.SetActive(true);
    }
}
