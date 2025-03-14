# Construct an image from the micromamba environment
FROM mambaorg/micromamba:2.0.5

COPY .. /app
WORKDIR /app

#RUN micromamba create --file /app/environment.yml -y
#RUN micromamba activate jsrat
