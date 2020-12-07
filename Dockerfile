FROM rocker/shiny-verse:3.6.3

RUN apt-get update && apt-get install -y libjpeg-dev

RUN echo "force3"
COPY r-requirements.dcf /r-requirements.dcf
RUN R -e 'source("https://gist.githubusercontent.com/cannin/6b8c68e7db19c4902459/raw/installPackages.R"); installPackages("/r-requirements.dcf")'

RUN cp -R /usr/local/lib/R/site-library/tumorcomparer/shinyApp/ /srv/shiny-server/tumorcomparer/

# Allows plotly to render
RUN chown -R shiny:shiny /srv/shiny-server

# Disables certain shiny-server protocols that prevent the app from loading at Dana-Farber and MD Anderson
# COPY inst/scripts/shiny-server.conf /etc/shiny-server/shiny-server.conf

CMD ["shiny-server"]
