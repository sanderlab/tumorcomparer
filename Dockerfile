FROM cannin/r-shiny-server:ubuntu-14.04.4_r-3.3.2_java-8_shiny-server-1.5.2.837

COPY inst/scripts/installPackage.R installPackage.R
RUN R -e 'source("installPackage.R")'

RUN cp -R /usr/local/lib/R/site-library/tumorcomparer/shinyApp/ /srv/shiny-server/tumorcomparer/

# Copy data from data folder
#COPY data/* /usr/local/lib/R/site-library/tumorcomparer/shinyApp/www/db/

# Allows plotly to render
RUN chown -R shiny:shiny /srv/shiny-server

# Disables certain shiny-server protocols that prevent the app from loading at Dana-Farber and MD Anderson
COPY inst/scripts/shiny-server.conf /etc/shiny-server/shiny-server.conf

CMD ["shiny-server"]
