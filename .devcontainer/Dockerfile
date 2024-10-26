FROM gcc:15

ARG USER=devcontainer
ARG PASSWD=""
ARG UID=1000
ARG GID=1000
ARG SHELL=/usr/bin/zsh

ARG DEBIAN_APT_REPOSITORY_URL="http://deb.debian.org"
# ARG DEBIAN_APT_REPOSITORY_URL="https://mirrors.tuna.tsinghua.edu.cn"
ARG LLVM_APT_REPOSITORY_URL="https://apt.llvm.org"
# ARG LLVM_APT_REPOSITORY_URL="https://mirrors.tuna.tsinghua.edu.cn/llvm-apt"

ARG GCC_MAJOR_VERSION="15"
ARG CLANG_MAJOR_VERSION="21"
ARG REINSTALL_CMAKE_VERSION_FROM_SOURCE="4.2.1"

RUN sed -i "s|http://deb.debian.org|${DEBIAN_APT_REPOSITORY_URL}|g" /etc/apt/sources.list.d/debian.sources
RUN sed -i "/^Components:/ s/$/ contrib non-free non-free-firmware/" /etc/apt/sources.list.d/debian.sources

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get -y install --no-install-recommends lsb-release python-dev-is-python3 \
    pinentry-curses netcat-openbsd ninja-build cppcheck gdb ccache valgrind autoconf-archive rsync \
    openssh-client gpg-agent locales tzdata curl wget git git-lfs sudo vim zsh xauth gnuplot

# Install specific version of cmake.
COPY ./reinstall-cmake.sh /tmp/
RUN if [ "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}" != "none" ]; then \
    chmod +x /tmp/reinstall-cmake.sh && /tmp/reinstall-cmake.sh ${REINSTALL_CMAKE_VERSION_FROM_SOURCE}; \
    fi
RUN rm -f /tmp/reinstall-cmake.sh

# Add LLVM-APT ppa and install LLVM compilers and libraries.
RUN bash -c "$(curl -fsSL ${LLVM_APT_REPOSITORY_URL}/llvm.sh)" "" ${CLANG_MAJOR_VERSION} all -m ${LLVM_APT_REPOSITORY_URL}

# Install flang compiler and library.
RUN apt-get -y install --no-install-recommends flang-${CLANG_MAJOR_VERSION} libflang-${CLANG_MAJOR_VERSION}-dev

# Clean apt cache
RUN apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Set the PATH for LLVM compilers.
RUN sed -i "$ a\export PATH=\"\/usr\/lib\/llvm-${CLANG_MAJOR_VERSION}\/bin\:\$\{PATH\}\"\nexport TERM=\"xterm-256color\"" /etc/bash.bashrc
RUN sed -i "$ a\export PATH=\"\/usr\/lib\/llvm-${CLANG_MAJOR_VERSION}\/bin\:\$\{PATH\}\"\nexport TERM=\"xterm-256color\"" /etc/zsh/zshrc

RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen
# RUN sed -i '/zh_CN.UTF-8/s/^# //g' /etc/locale.gen
RUN locale-gen
RUN update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

RUN groupadd -g ${GID} ${USER}
RUN useradd -s ${SHELL} -u ${UID} -g ${GID} -G sudo -m ${USER}

RUN if [ -n "${PASSWD+x}" ] && [ -n "$PASSWD" ]; then \
    echo "root:$PASSWD" | chpasswd; \
    echo "$USER:$PASSWD" | chpasswd; \
    else \
    passwd -d root; \
    passwd -d "$USER"; \
    fi
