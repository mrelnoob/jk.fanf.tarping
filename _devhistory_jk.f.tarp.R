############## *-------------------------------------------------* #############
####################### DEV_HISTORY FOR {jk.fanf.tarping} ######################
############## *-------------------------------------------------* #############

# All comments in this file or in this project are first meant for future me and, second,
# for other so they can reproduce my work. That is why there are so many of them, so anyone can
# do what I did, regardless of his/her own experience. Note also that English and R are not my
# native tongues, so please forgive any mistakes.
# IMPORTANT NOTE: I chose to use a {targets} pipeline to organise my data analysis procedure to ensure
# a clean workflow and favour reproducibility. However, that also means that my code cannot be run as is
# if you do not know how to use {targets}, so be careful!
# Note also that I briefly explain how to use {targets} at the end of this document and I will add
# {renv} later, if I can.





# -------------------------------- #
##### 0. Important R reminders #####
# -------------------------------- #

##### * 0.1. NA handling -------------------------------------------------------
# ---------------------------------------------------------------------------- #

# Although it may seem strange, R doesn't handles NA as we could expect it. The use of "is.na()" or
# "!is.na()" should be preferred over the use of classical equality operators "==" or "!=".
# Cf. https://stackoverflow.com/questions/28857653/removing-na-observations-with-dplyrfilter



##### * 0.2. Factor and character variables ------------------------------------
# ---------------------------------------------------------------------------- #

# R is kind of a pain in the arse when it comes to dealing with categorical variables. Theoretically,
# these are `character` objects but most R functions expect `factor` objects. So you have to frequently
# juggle between those two types of data format and remember that a `factor` has `levels` and thus that
# what you see on your screen is not how R handles categorical variables internally. So to modify a
# `factor` object, you often have to convert it first into a `character` object first.
# Cf. many of my custom functions.



##### * 0.3. Package building, {targets} and RMarkdown -------------------------
# ---------------------------------------------------------------------------- #

# To optimally associate Package building, {targets}, renv and RMarkdown, I should probably refer to:
# https://books.ropensci.org/targets/





########################## ********************************************* ###############################
# ------------------------- #
##### 1. Project set up #####
# ------------------------- #

##### * 1.1. Keeping track of changes (Git) ------------------------------------
# ---------------------------------------------------------------------------- #

# To keep track of all future changes and have a backup, I need to initiate a Git version control
# repository and link my R project folder to my GitHub account. Here, I created the Git repository
# alongside the R project (using the buttons of RStudio), but I could have done it with the procedure
# described below. As it is not the first time I do that on my computer and with this R version, my Git
# is already parametrised and I do not need to follow the entire procedure.
# As I'm connecting my 'jk.fanf.tarping" project AFTER its creation, I'm following a 'GitHub-last' procedure:
# STEP 1: Make sure Git is installed (e.g. by typing 'git status' in the R terminal), if not, install it!
# STEP 2: Enter your project and type 'git init' in the R Terminal, it will initiate a Git
# repository within your project:
system("git init") # NOTE: the 'system()' function enables sending commands to the Terminal/CLI,
# that is to work in "command-lines" (hence the name CLI, for "command-lines interface").
# It works yet the problem is that RStudio seems unable to display the CLI responses like a
# regular CLI, preventing interactions with it. That is why I sometimes say that I had to do
# stuff directly in a CLI/Terminal.
# STEP 3: Verify your configuration:
system("git config --global user.name 'fanf'")
system("git config --global user.email 'francois.martin74190@gmail.com'")
# Of course, you need to personalize the user.name and email. If you want the config to be
# only true for the current project (and not global), remove the '--global' from the previous
# lines.
# STEP 4: If you've never done it before (on this computer), you also need to set-up SSH keys
# to be able to connect to Git and GitHub without supplying you username and password every
# time you do something.
# The full procedure is described here: https://help.github.com/articles/generating-ssh-keys
# But here goes:
system("ls -al ~/.ssh") # Check for existing SSH keys. If none is found, you must generate one.
system("ssh-keygen -t rsa -C 'francois.martin74190@gmail.com'") # Generates a new SSH key
# using the rsa algorithm (however, I was forced to do that directly in a CLI).
system("eval '$(ssh-agent -s)'") # Checks if the ssh-agent is launched (same --> CLI).
# If it works, it should return something like "Agent pid 59552".
system("ssh-add ~/.ssh/id_rsa") # Adds the new SSH private key to the ssh-agent! (CLI).
# Then you have to add the public key to your GitHub account settings:
# e.g. https://kbroman.org/github_tutorial/pages/first_time.html
system("ssh -T git@github.com") # To check if it works. If it does, it should answer with
# something like "Hi mrelnoob! You've successfully authenticated, but Github does not provide
# shell access". There is clearly a bug in RStudio (or R) that prevents me from doing all
# that from RStudio, but anyway...

# STEP 5: Make some changes in the project and make your first commit:
usethis::use_git(message = ":tada: Initial commit") # Then restart RStudio and the Git tab will
# appear on the top-right box of RStudio.
# STEP 6: Log in your GitHub account and create a new repository (without anything in it).
# STEP 7: Use the following command to associate your R project with the GitHub project:
system2("git remote add origin git@github.com:mrelnoob/ppl.ipa") # Here also, personalize
# with your own account and project names! And here again, it does not work (so --> CLI).
# STEP 8: Finally, you can push the changes made to your local files to GitHub:
system2("git push -u origin master") # Same (CLI).
# Even using a CLI (e.g. GitBash), you may receive an error message saying that the remote
# repository contains work that you do not have locally (...). It happens if you included files
# in your GitHub project when you created it (e.g. a README or a LICENCE file). So theoretically,
# you need to always start with a "pull" before you push anything! If, like me, you forgot,
# you'll be in trouble and won't be able to pull. To force Git, you may use "git push -f origin
# master" (the -f means to "force" the push).
# IMPORTANT NOTE: because of the CLI-RStudio bugs, I can "commit" from RStudio but I cannot push,
# so I will always be forced to do it from a CLI every time!

# To ignore changes made to the Rproj file:
usethis::use_git_ignore("jk.fanf.tarping.Rproj")



##### * 1.2. Managing the project library (renv) -------------------------------
# ---------------------------------------------------------------------------- #
renv::init() # To initiate the project local library manager named 'renv'. Among other things, it
# created a 'renv' folder and a 'renv.lock' file. This file is the only thing that you need to be
# able to recreate the same R environment (correct package versions) as the one I used in this
# project PROVIDED THAT you use the same R version as me (currently R 4.2.2).

# To install a package into the 'renv' project, use:
renv::install() # You can either specify the package name manually (e.g. "ggplot2") or you can call
# the function as is and it will install all the packages of your DESCRIPTION file (if you have one),
# and also screen R and Rmd files and install all packages mentioned as 'library(pkg)', 'require(pkg)',
# and 'pkg::fun()'! Pretty cool, right?

# To check the status of your 'renv' project:
renv::status()

# To update the lockfile (to add or remove packages):
renv::snapshot() # It will compare what's inside the R project and what is listed in the lockfile and
# ask if you want to add or remove the packages that are not yet on both projects. You can also use:
renv::clean()

# To restore another 'renv' project (if you want to create the local library of another project, yours
# or someone else's):
renv::restore() # You may need to use 'renv::init()' first if you only have a lockfile. Remember though,
# you need to install the same major R version as the one used to create the 'renv' project!

# Other functions:
renv::deactivate()       # Deactivate local environment
renv::activate()         # (Re)activate local environment
renv::dependencies()     # List used packages (R and Rmd files)
renv::history()          # Browse previous commits (with git)



##### * 1.3. Project architecture ----------------------------------------------
# ---------------------------------------------------------------------------- #

# To create a folder containing my data and functions:
dir.create("data")
dir.create("R")
# To create other useful folders for my project:
dir.create("output")
dir.create("output/plots")
dir.create("output/tables")
dir.create("output/texts") # I could also create my packages directly in my functions, that would be
# particularly appropriate if I implement a {target} pipeline.

# If I am building an R package, then these folders are not supposed to be here and should be ignored:
# usethis::use_build_ignore("output/")
# usethis::use_build_ignore("tables/")
# usethis::use_build_ignore("texts/")
# usethis::use_build_ignore("plots/")
# usethis::use_git_ignore("plots/") # To avoid saturating Git, I ignore the folders prone to contain
# # rather heavy files such as spatial layers and plots, but not tables and texts!
# # NOTE: the best way to build a package is not this one (see other '_devhistory' files or refer to
# # section 0.3. of this document).



##### * 1.3. Creating scripts for custom functions -----------------------------
# ---------------------------------------------------------------------------- #

file.create(... = "R/01_01_prep_functions.R")
file.create(... = "R/01_02_cleaning_functions.R")
file.create(... = "R/01_03_exploration_functions.R")








##### A CONTINUER§§§§§ ----
##### A CONTINUER§§§§§ ----
##### A CONTINUER§§§§§ ----

##### * 1.4. Creating reports (RMarkdown) --------------------------------------
# ---------------------------------------------------------------------------- #

file.create(... = "output/texts/ppl.tits.exploration_report.Rmd")
file.create(... = "output/texts/ppl.tits.intermediate_analyses_report.Rmd") # Using this command,
# a .Rmd file will be created but will lack the YAML header skeleton that should thus be manually
# placed at the top of the document.
# NOTE: for guidance to use .Rmd documents with {targets}, please refer to the {targets} related
# section (chapter 2 below).


### ** 1.3.1. To manage citations and bibliography ----
# _____________________________________________________

# To manage citations and get an automatic bibliography with RMarkdown, I have to follow these
# steps:
#  1) Using Zotero (or something similar), I have to 'export' the references to be cited in the
#     report in a BibTex format (.bib) and place this text file in the same folder as my .Rmd file.
#  2) Call this document in the `bibliography` field in the YAML metadata header (e.g.
#     bibliography: my_example.bib).
#  3) In text, I use arobases (@) and brackets ([], use semi-colons ";" for separation between
#     references) to add citations (e.g. "@Martin2022 said that..." or "blabla [@Martin2022;
#     see also @Darwin1832]").
#  4) I can change the citation style by using the `csl` field in the YAML metadata header
#     (e.g. csl: my_style.csl) and pasting the said style in the same folder as before.
# Thus, I pasted my BibTex file in the same folder as my .Rdm file. But in the case where my .Rmd
# file would be at the root of my package, I need to tell R to ignore it:
usethis::use_build_ignore("ppl.tits_biblio.bib")
# For practical reasons, this .bib file will certainly be updated many times during the duration
# of the project. Also, it may be useful to manually edit the file to shorten the reference tags
# since Zotero tends to create long tag using the name of the 1st author, the 1st work of the title
# and the year of publication.





########################## ********************************************* ###############################
# -------------------------------- #
##### 2. Populating my project #####
# -------------------------------- #

##### * 2.1. Writing functions -------------------------------------------------
# ---------------------------------------------------------------------------- #

# I may now write my scripts and associated functions --> THEN TARGETS§§§§§
# I may now write my scripts and associated functions --> THEN TARGETS§§§§§
# I may now write my scripts and associated functions --> THEN TARGETS§§§§§
# I may now write my scripts and associated functions --> THEN TARGETS§§§§§




########### *-----------------------------------------------------* ############
############################ Main Git commits ##################################
# ---------------------------------------------------------------------------- #
usethis::use_git(message = ":boom: Created the project structure!")
usethis::use_git(message = ":metal: Created a new function")
usethis::use_git(message = ":zap: Ignoring something")
usethis::use_git(message = ":pencil: Documented a function or wrote something")
usethis::use_git(message = ":hammer: Ongoing programming!")
usethis::use_git(message = ":white_check_mark: Updated the {target} pipeline")
usethis::use_git(message = ":x: Problem detected!")
#system("git push") # Or using a CLI!
# Don't forget to push your commits once you're sure you made no mistakes.
# ---------------------------------------------------------------------------- #
# ------------------------------- THE END ------------------------------------ #
########### *-----------------------------------------------------* ############
