     sudo apt-get update --fix-missing
     sudo rm -rf /var/lib/apt/lists/*
     sudo rm -rf /var/lib/apt/lists/partial/*
     sudo apt-get clean
     sudo apt-get update
     sudo apt-get update -o Acquire::No-Cache=True
     sudo apt-get update -o Acquire::BrokenProxe=True
     sudo rm -vf /var/lib/apt/lists/*     
     sudo apt-get update     -y                                    
     sudo apt-get upgrade   -y
     sudo apt-get dist-upgrade  -y 
     sudo apt autoremove -y  

