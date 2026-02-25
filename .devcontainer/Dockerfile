ARG BASE_IMAGE="gcc:latest"
FROM ${BASE_IMAGE}

ARG GCC_MAJOR_VERSION
ARG CLANG_MAJOR_VERSION
ARG CMAKE_VERSION

ARG USER=devcontainer
ARG PASSWD=""
ARG UID=1000
ARG GID=1000
ARG SHELL=/usr/bin/zsh
ARG TZ="UTC"

ARG DEBIAN_APT_REPOSITORY_URL="http://deb.debian.org"
# ARG DEBIAN_APT_REPOSITORY_URL="https://mirrors.tuna.tsinghua.edu.cn"
ARG LLVM_APT_REPOSITORY_URL="https://apt.llvm.org"
# ARG LLVM_APT_REPOSITORY_URL="https://mirrors.tuna.tsinghua.edu.cn/llvm-apt"

RUN export DEBIAN_FRONTEND="noninteractive" \
    && sed -i "s|http://deb.debian.org|${DEBIAN_APT_REPOSITORY_URL}|g" /etc/apt/sources.list.d/debian.sources \
    && sed -i "/^Components:/ s/$/ contrib non-free non-free-firmware/" /etc/apt/sources.list.d/debian.sources

RUN apt-get update && apt-get -y install --no-install-recommends lsb-release python-dev-is-python3 \
    pinentry-curses netcat-openbsd ninja-build cppcheck gdb ccache valgrind autoconf-archive rsync \
    openssh-client gpg-agent dialog locales tzdata curl wget zip unzip tar git git-lfs sudo vim zsh \
    xauth gnuplot \
    && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/* \
    && sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && sed -i '/zh_CN.UTF-8/s/^# //g' /etc/locale.gen \
    && locale-gen && update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 \
    && sed -i "$ a\export TERM=\"xterm-256color\"" /etc/bash.bashrc \
    && sed -i "$ a\export TERM=\"xterm-256color\"" /etc/zsh/zshrc

# Install specific version of cmake.
RUN if [ -n "${CMAKE_VERSION}" ]; then \
    apt-get update \
    && bash -c "$(curl -fsSL https://raw.githubusercontent.com/devcontainers/templates/main/src/cpp/.devcontainer/reinstall-cmake.sh)" "" ${CMAKE_VERSION} \
    && rm -rf /var/lib/apt/lists/*; \
    fi

# Add LLVM-APT ppa and install LLVM compilers and libraries.
RUN apt-get update && bash -c "$(curl -fsSL ${LLVM_APT_REPOSITORY_URL}/llvm.sh)" "" ${CLANG_MAJOR_VERSION} all -m ${LLVM_APT_REPOSITORY_URL} \
    && apt-get -y install --no-install-recommends flang-${CLANG_MAJOR_VERSION} libflang-${CLANG_MAJOR_VERSION}-dev \
    && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/* \
    && sed -i "$ a\export PATH=\"\/usr\/lib\/llvm-${CLANG_MAJOR_VERSION}\/bin\:\$\{PATH\}\"" /etc/bash.bashrc \
    && sed -i "$ a\export PATH=\"\/usr\/lib\/llvm-${CLANG_MAJOR_VERSION}\/bin\:\$\{PATH\}\"" /etc/zsh/zshrc

RUN groupadd -g ${GID} ${USER} && useradd -s ${SHELL} -u ${UID} -g ${GID} -G sudo -m ${USER} \
    && if [ -n "${PASSWD}" ]; then \
    echo "root:${PASSWD}" | chpasswd; \
    echo "${USER}:${PASSWD}" | chpasswd; \
    else \
    passwd -d root; \
    passwd -d "${USER}"; \
    fi

USER ${USER}
RUN zsh -c "$(curl -fsSL https://raw.githubusercontent.com/BoyleDevTeam/boyle/master/.devcontainer/configure-user.zsh)"
