FROM python:3.10-slim-bullseye

WORKDIR /home/app

RUN pip install ngstrefftz 
ENV PYTHONPATH=$PYTHONPATH:/usr/local/lib/python3.10/site-packages/
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

RUN pip3 install numpy webgui_jupyter_widgets notebook
#RUN pip3 install jupyter_contrib_nbextensions 
#RUN jupyter nbextension enable --py widgetsnbextension
#RUN jupyter nbextension enable --py webgui_jupyter_widgets

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}
#ENV mv /home/app/ngstrefftz/* /home/${NB_USER}/

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

USER ${NB_USER}

WORKDIR /home/${NB_USER}
#RUN git clone https://github.com/PaulSt/NGSTrefftz ${HOME}/ngstrefftz
COPY ./docs/notebooks ${HOME}/docs/notebooks

CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root", "docs/notebooks" ]  
