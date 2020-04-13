# Based on:
* http://blog.gopheracademy.com/advent-2014/easy-deployment/
* https://github.com/ArturoBlas/docker-captainhook

# Pull captainhook Docker image
sudo docker pull cannin/captainhook

# Run captainhook Docker image
sudo docker run -d -p 8080:8080 -v /webhooks:/webhooks cannin/captainhook

# Run from local installation
* https://github.com/bketelsen/captainhook
captainhook -listen-addr=0.0.0.0:3830 -echo -configdir /webhooks &

# Trigger captainhook
* NOTE: Must run: chmod +x kidneyMetabProject.sh
curl -X POST http://127.0.0.1:3830/kidneyMetabProject
