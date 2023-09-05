FROM tomcat:9.0
COPY ./JalTantra-Website/classes/artifacts/Jaltantra/Jaltantra.war /usr/local/tomcat/webapps/ 
EXPOSE 8080
