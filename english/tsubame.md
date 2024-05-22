# TSUBAME
* In the following, we assume as the user is the Tokyo-tech student or employee.

# Making account
1. Login to the Tokyo-tech portal, and click "TSUBAME portal"
2. Fill the required fields
3. E-mail comes in moment, and the click the URL in the mail
4. Confirm the account name

# ssh key generation
1. Open the terminal in mac/windows-wsl and type `ssh-keygen`
2. Specify the name of public key
3. Set your passphrase if you like
4. Copy the texts in the above public key file
5. Back to the TSUBAME portal, and paste the copied text into the field of "SSH public key registration"
6. Click "add"
7. Go back to terminal, and you can login like `ssh [account_name]@login.t3.gsic.titech.ac.jp -i [private_key-file]`
8. Setup the ssh configure file `~/.ssh/config` to make login easier. An example of minimal configuration is,
    ```bash
    Host tsubame
        HostName login.t3.gsic.titech.ac.jp
        User your_name
        IdentityFile /Users/your_name/.ssh/id_rsa_tsubame
    ```
9. You can login TSUBAME by `ssh tsubame`.

# Getting computational resources
* You need to make a group to get computational resources.
* This section involves the budget payment procedures, so often only instructors need to manage them.

## Making group
1. login to "TSUBAME-portal"
2. click the making group, and fill the required fields
3. you can manage the groups with "managing the group" -> "detail"

## Payment procedure
### Registering purchase code
1. you need the "purchase code" for paying the points
2. click "managing the purchase code" --> "apply for the new purchase code"
3. fill the required fields
4. wait
5. the purchase code will be issued, and then fill the budget information from "managing the group"

### Purchasing points
1. login to the TSUBAME portal
2. click "managing the group"
3. select the belonging groups
4. click "detail"
5. purchasing points

# Submitting jobs
## Normal
* You can execute your calculation as *jobs* to the supercomputer.
* Supercomputer queuing system takes care of your (and also others') jobs.
* To register your job, execute: `qsub -g [group_name] script.sh`
    * If you don't specify the group name, the job will be a trial-run so it stops in 1 hour.
* The `script.sh` file contains the procedure of your calculation.
* An example of `script.sh` for running Gaussian is (input file: `h2o.com`)
    ```bash
    #!/bin/bash
    #$ -cwd
    #$ -l f_node=1
    #$ -l h_rt=00:10:00
    #$ -V

    . /etc/profile.d/modules.sh
    module load gaussian16

    g16 h2o.com
    ```

<!--
## Using booked node
* You can book the nodes via TSUBAME portal.
* With booked nodes (AR_ID should be given), you can submit jobs by
`qsub -g [group_name] -ar [AR_ID] script.sh`.
-->

# Confirming jobs
* To confirm your job status, type: `qstat`.
* Job states
    * `qw`: waiting for run
    * `r` : running
    
# Stopping jobs
* To stop your jobs, type: `qdel [JOB_ID]`.
* To know the JOB_ID, do `qstat`.