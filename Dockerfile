FROM paulstdocker/ngstrefftz:latest

##RUN pip3 install jupyter
#RUN pip3 install jupyter_contrib_nbextensions
#RUN pip3 install jupyter_nbextensions_configurator
#RUN pip3 install RISE
#RUN pip3 install ipywidgets
#RUN jupyter contrib nbextension install 
##RUN cd /home/app/ngstents/tentswebgui && pip3 install --user .

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

USER ${NB_USER}


RUN jupyter nbextensions_configurator enable --user
RUN jupyter nbextension enable --user codefolding/main 
RUN jupyter nbextension enable --user --py tentswebgui

#WORKDIR /home/${NB_USER}
#COPY . ${HOME}
WORKDIR /home/app/ngstrefftz/docs/notebooks

CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root" ]  
