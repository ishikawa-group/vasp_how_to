# TSUBAME
* In the following, we assume as the user is the Tokyo-tech student or employee.

## Making account
1. login to the Tokyo-tech portal, and click "TSUBAME portal"
2. fill the required fields
3. e-mail comes in moment, and the click the URL in the mail
4. confirm the account name, and registar the ssh key
    4-1. open the terminal in mac/windows-wsl and type `ssh-keygen`
    4-2. specify the name of pubilc key
    4-3. copy the texts in the above public key file
    4-4. back to the TSUBAME portal, and paste the copied text into the field of "SSH public key registration"
    4-5. click "add"
5. go back to terminal, and you can login like `ssh [your_account_name]@login.t3.gsic.titech.ac.jp -i [your_private_keyfile]`
6. setup the ssh configure file to make login easier

## Getting resources
* You need to make a group to get computational resources.

### Making group
1. login to "TSUBAME-portal"
2. click the making group, and fill the required fields
3. you can manage the groups with "managing the group" -> "detail"

### Payment procedure
#### Registring purchase code
1. you need the "purchase code" for paying the points
2. click "managing the purchase code" --> "apply for the new purchase code"
3. fill the required fields
4. wait
5. the purchase code will be issued, and then fill the budget information from "managing the group"

#### Purchasing points
1. login to the TSUBAME portal
2. click "managing the group"
3. select the beloging groups
4. click "detail"
5. purchasing points

##

## Submitting jobs
* `qsub -g [group_name] script.sh`
* If you don't specify the group name, the job will be a trial-run.

## Confirming jobs
* `qstat`