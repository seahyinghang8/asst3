# Google Cloud Platform Setup Instructions #

In order to test your assignment for performance points, you will need to run it on a VM instance on the Google Cloud Platform (GCP). We've already sent you student coupons that you can use for billing purposes. All that is left is for you to create your own VM in the cloud. Here are the steps for how to get setup for running on GCP.

NOTE: For those working in teams, it might be desirable for both students to use the same virtual machine. To do so, only one of you should first create the VM instance following the instructions below. Then, follow the instructions at https://cloud.google.com/compute/docs/access/granting-access-to-resources to grant access to your partner.

NOTE #2: __Please don't forget to SHUT DOWN your instances when you're done with your work for the day!  The VM Instances you create cost approximately $1.13 per hour.  Leaving it on accidentally for a day could quickly eat up your $50 per student quota for the assignment.__
 
  1. Go to https://google.secure.force.com/GCPEDU?cid=GOWuE1bzxv%2BQeeqn4jG9xArix82DMPXlZ%2FuplPh3CAY6Uye5YnluvOqYiCjFizqY, and fill out the form with your personal info.
  
  2. Once you receive a coupon code, log into your Gmail account, go to https://console.cloud.google.com/education, and put in the code. Click accept. __Note:__ If you're using your Stanford account, you will probably get the message `"You don’t have permission to create a billing account for your organization."` In that case, you should use a personal Gmail account instead (e.g. `myemail@gmail.com`).
     
  3. Now, go back to the dashboard by clicking the `Google Cloud Platform` tab. Create a new project by first clicking on the upper tab `Select a new project` and then clicking `new project` on the newly opened window. Make sure to select your newly-created project now.
 
  4. By default, GCP limits new users to using machines with at most 24 vCPUs, but your assignment requires you to run performance tests on 32 vCPUs. (Note: a GCP vCPU is equivalent to a hyper-thread.)  Before proceeding, you will need to request a quota increase to at least 32 vCPUs. Click on the Navigation Menu button (located on the upper-left side of the page), and under `IAM & Admin`, select `Quotas`. Click `Edit Quotas` and click the box by `Compute Engine API (us-west-1)`.  On the right hand side of the screen you will see menus that ask you to enter your name and email, and can request a quota increase.  Please request up to 32 vCPUs.  GCP wil almost certainly automatically grant your request.

  5. Click on the Navigation Menu button (located on the upper-left side of the page) and you will get a list of services. 
  
  6. Click on `Compute Engine`, and then click on `VM Instances`. 
  
  7. Click on `enable billing`. This will allow you to use your Google Cloud Patform Credits. 
  
  8. Now you're ready to create a VM instance. Click on the button that says `Create Instance`. Fill out the form such that your cloud-based VM has the following properties: 
       - Region __us-west1__ (Oregon)
       - Type n1-highcpu-32 (__32 vCPUs__, __28.8 GB__ memory) 
       - Ubuntu __18.04 LTS__  
       - At least a __20GB__ Standard persistent disk.

  9. Finally, click `Create` 
  
Now that you've created a VM instance on Google Cloud, you should be able to __SSH__ into it. 2 ways to do so are to either open it from the VM instance page (the one that shows a list of all your VM instances) or use the `gcloud` command (which you can use it on your local console or on the Google Cloud shell). The `gcloud` command would look like this: 

      gcloud compute instance MY_INSTANCE_NAME

You can find `MY_INSTANCE_NAME` in the *VM instances* page.

Once you finally SSH into your VM instance, follow the instructions below to finish setting up your cloud environment:

  1. Install the following programs using the `sudo apt-get install` command:
      - g++ 
      - git 
      - make
      - your desired editor (emacs/vim)
      
  To be more specific, you should use these commands in order:  
     
        sudo apt-get update
        sudo apt-get install g++ git make emacs
        
  2. Git clone the assignment repo 
  
  3. Download the graph files from the class directory using the following command: 
  
    scp YOURSUNETID@myth.stanford.edu:/afs/ir.stanford.edu/class/cs149/data/asst3_graphs/all_graphs.tgz .
  
  4. Build and run the starter code as indicated by the assignment handout. __Note: When running on GCP, you will need to change the variable `REF_LIB` in the `pagerank` and `bfs` Makefile's so that it points to the reference solution binary that was precompiled for GCP (and not myth).__ 
  
If you're confused about any of the steps, having problems with setting up your account or have any additional questions, reach us out on Piazza!
  
__Again, please don't forget to SHUT DOWN your instances when you're done with your work for the day!__
