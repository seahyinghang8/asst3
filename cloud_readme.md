# Cloud Setup Instructions #

In order to test your assignment, you will need to run it on a Google Cloud VM Instance. We've already sent you student coupons that you can use for billing purposes. All that's left is for you to create your own VM in the cloud! Here are the steps for how to do so: 

  1. Go to https://google.secure.force.com/GCPEDU?cid=GOWuE1bzxv%2BQeeqn4jG9xArix82DMPXlZ%2FuplPh3CAY6Uye5YnluvOqYiCjFizqY, and fill out the form with your personal info.
  
  2. Once you receive a coupon code, log into your Gmail account, go to https://console.cloud.google.com/education, and put in the code. Click accept. __Note:__ If you're using your Stanford account, you will probably get the message `"You don’t have permission to create a billing account for your organization."` In that case, you should use a personall Gmail account (e.g. `myemail@gmail.com`).
  
  3. Now, go back to the dashboard by clicking the `Google Cloud Platform` tab. Create a new project by first clicking on the upper tab `Select a new project` and then clicking `new project` on the newly opened window. Make sure to select your newly-created project now.
  
  4. Click on the Navigation Menu button (located on the upper-left side of the page) and you will get a list of openings. 
  
  5. Click on `Compute Engine`, and then click on `VM Instances`. 
  
  6. Click on `enable billing`. This will allow you to use your Google Cloud Patform Credits. 
  
  7. Now you're ready to create a VM instance. Click on the button that says `Create Instance`. Fill out the form such that your cloud-based VM has the following properties: 
       - Region us-west1 (Oregon)
       - Type n1-highcpu-32 (32 vCPUs, 28.8 GB memory) 
       - Ubuntu 18.04 LTS 
       - At least a 20GB Standard persistent disk.

  8. Finally, click `Create` 
  
Now that you've created a VM instance on Google Cloud, you should be able to __SSH__ into to. 2 ways to do so are to either use the `gcloud` command (which you can use it on your local console or on the Google Cloud shell) or to open it from the VM instance page (the one that shows a list of all your VM instances). 

Once you've finally SSH'd into your VM instance, you need to do the following so that it's ready to run your programming assignment: 

  1. Install the following programs using the `sudo apt-get installl` command:
      - g++ (version >= 5.0)
      - git 
      - make
      - your desired editor (emacs/vim)
      
  2. Git clone the assignment repo 
  
  3. Download the graph files from: `myth.stanford.edu:/afs/ir.stanford.edu/class/cs149/data/asst3_graphs/all_graphs.tgz`
  
  4. Build and run the starter code as indicated by the assignment handout. __Note:__ you need to change  REF_LIB in the MakeFile's so that it points to the cloud reference binary. 
  
  
__Please Don't Forget to SHUT DOWN your instances when you're done with your work for the day!__
