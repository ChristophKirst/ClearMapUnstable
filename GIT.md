git
===


Install
-------

plain user:

 * in terminal execute: 
 
    `cd basedirectory`
 
    `git clone https://git.assembla.com/idisco.git`


developer / programmer:

  * create account at assembla

  * got to [https://git.assembla.com/idisco.git](https://git.assembla.com/idisco.git) and press fork button 

  * in terminal execute:
	
	`cd basedirectory`

	`git clone https://git.assembla.com/idisco.git`
	
  * configure remotes (named upstream)
        
	`cd iDisco`

	`git remote add upstream https://git.assembla.com/idisco.git`

	`git fetch upstream`


Backup
------

to backup your version in case you followed the developer / proogrammer route:

  * in terminal in the iDisco directory execute:

      `git add -A`

      `git commit -m 'some description of what you did'`

      `git push`


Update
------    

plain user:

  * in terminal in the iDisco directory execute
     
      `git pull`


programmer: 

in case you want to update your code from the upstream repository

  * in terminal execute:
 
      `git fetch upstream`
      
      `git merge upstream/master`

  * if mergin fails, some files will be highlighted with <<<<<< >>>>>> entries, fix this manually

  * if you dont care about your own changes and simply want the plain new version:

      `git reset --hard upstream/master`

  * to force it to your fork on github use
       
	  `git push origin/master --force` 


Submitting
----------

in case you have something to contribute to the code:
 
  * follow the steps in the Backup section first

  * got to [https://git.assembla.com/idisco.git](https://git.assembla.com/idisco.git and click pull request 


Refs
----

iDisco home:

  * [https://github.com/ChristophKirst/iDisco](https://github.com/ChristophKirst/iDisco)

A good source to get questions answered: 

  * [https://help.github.com/](https://help.github.com)
  * [http://git-scm.com/documentation](http://git-scm.com/documentation)

git home:

  * [http://git-scm.com](http://git-scm.com/)

