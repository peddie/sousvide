{
  description = "Development shell for sous vide modelling calculations.";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils/master";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        rPackages = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            # Tidyverse
            tidyverse
            modelr
            data_table
            reshape2
            slider

            # Stan stuff
            rstan
            brms
            pkgs.cmdstan

            # Bayesian stuff
            tidybayes
            bayesplot
            dagitty
            coda
            posterior

            # Visualization
            ggplot2
            raster
            cowplot
            gridExtra
            ggrepel
            RColorBrewer
            leaflet
            hexbin

            # Linear algebra
            RSpectra
            EigenR
            RcppEigen
            Rcpp
            Rcpp11
            dqrng

            # Spatial stuff
            sp
            sf
            rgdal

            # Disk I/O
            qs

            # Output
            knitr
            roxygen2
            docstring

            # Development
            lintr
            reprex
            roger
            styler
            tictoc
            lobstr
            microbenchmark
            devtools
            # TODO(MP): This package ends up with C++ linker errors.
            (buildRPackage {
              name = "cmdstanr";
              src = pkgs.fetchurl {
                url = "https://mc-stan.org/r-packages/src/contrib/cmdstanr_0.5.3.tar.gz";
                sha256 = "DF3BauB8DT7M9AEth3BXZlF0zuxosnKYxF+HPhRmbCM=";
              };
              propagatedBuildInputs = [ checkmate data_table jsonlite posterior processx R6 ];
            })
          ];
        };
        packageName = "sousvide";
      in {
        defaultPackage = pkgs.R;

        devShell = pkgs.mkShell {
          buildInputs = with pkgs; [
            git
            glibcLocales
            openssl
            which
            openssh
            curl
            wget
            file
            rPackages
          ];
          shellHook = ''
    mkdir -p "$(pwd)/_libs"
    export R_LIBS_USER="$(pwd)/_libs"
    echo ${rPackages}/bin/R
  '';
        };
      });

}
