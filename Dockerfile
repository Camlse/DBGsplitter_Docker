FROM Debian 3.16.51-3
COPY ./install.sh ./
rum chmod +x ./install.sh && ./install.sh
