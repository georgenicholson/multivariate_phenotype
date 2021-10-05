# Can configure username
git config --global user.name "georgenicholson"

# # Clone files from origin (not sure how often need to do this)
# git clone https://github.com/georgenicholson/multivariate_phenotype.git

# cd to Git controlled directory
cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/github_multivariate_phenotype/multivariate_phenotype
# Look at local branches
git branch
# switch to local main branch
git checkout main
# check everything up to date
git remote update && git status 
# Make the local master sync with remote master
git pull origin main

git merge main

# Checkout a local branch to work on
git checkout -b george
# Switch to local george branch
git checkout george
# git checkout main

#################################
# Make changes here
#################################

# Switch to local george branch
git checkout george
git add -A # Add all files to be committed
git commit -m 'Initial organisation 2' # Commit files with message
git push --set-upstream origin george

# switch to local main branch
git checkout main
# check everything up to date
git status
# delete local george branch
git branch -d george  
# Look at local branches
git branch

# delete branch at origin
# git push origin --delete george



# git diff george main

# Command for helping with credentials to avoid typing password each time
# git config --global credential.helper store

# git diff origin/main
# git diff HEAD origin/main

# 
# git reset HEAD *.RData
# git reset HEAD **/.Rproj.user/**

# git config --global http.proxy http://username:password@proxy.server.com:8080
# git config --global https.proxy http://username:password@proxy.server.com:8080
# 
# //Replace username with your proxy username
# //Replace password with your proxy password
# //Replace proxy.server.com with the proxy domain URL.
# //Replace 8080 with the proxy port no configured on the proxy server.
# 
# 
# 31882a5c1285c201d07778815fe9156599da419d 
# 
# git help -a | grep credential-store
# git help credential-store
# 
# git remote -v
# 
# git clone https://github.com/georgenicholson/covid_sampling
# 
# cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_sampling_github/covid_sampling
# git remote -v
# # Switching remote URLs from HTTPS to SSH
# git remote set-url origin git@github.com:georgenicholson/covid_sampling
# git remote -v
# eval "$(ssh-agent -s)"
# ssh-keygen -t ed25519 -C "nicholso@stats.ox.ac.uk"
# ssh-add ~/.ssh/id_ed25519
# clip < ~/.ssh/id_ed25519.pub
# ssh-add
# ssh-add -l -E sha256
# ssh -vT git@github.com
# cd /mnt/c/Users/nicho/Documents/bauer_sync/projects/covid/covid_sampling_github/covid_sampling
# # git config credential.helper store
# git push git@github.com:georgenicholson/covid_sampling origin main
# git add .
# git status -s
# git commit -m "Added folders"
# git push origin main
# git rev-parse --show-toplevel
# 
# 
# git push 'https://georgenicholson:github007@github.com/georgenicholson/covid_sampling'  origin main
# sudo apt install geomview
# 
# 
