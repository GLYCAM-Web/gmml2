FROM debian:latest

# For tzdata
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
#tzdata && \
#ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime && \
#dpkg-reconfigure --frontend noninteractive tzdata && \
#apt-get install -y \
git \
cmake \
g++ \
libgomp1 \
&& rm -rf /var/lib/apt/lists/*;

# Set working directory
WORKDIR /app/

# Clone the GitHub repository
RUN git clone https://github.com/GLYCAM-Web/gmml2.git gmml2
WORKDIR /app/gmml2/

#COPY ./make.sh ./compilePrograms.sh ./CMakeLists.txt ./version.txt ./versionHeader.cmake /app/gmml2 
#COPY ./cmakeFileLists /app/gmml2/cmakeFileLists
#COPY ./src /app/gmml2/src 
#COPY ./include /app/gmml2/include 
#COPY ./programs /app/gmml2/programs
#COPY ./bin/ /app/gmml2/bin
#COPY ./tests/ /app/gmml2/tests
RUN ./make.sh -j10




# Optional but allows you to move exe files.
#ENV LD_LIBRARY_PATH=/app/gmml2/lib/:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH=/app/gmml2/lib

