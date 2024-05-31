{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
  };

  outputs = inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      # This is the list of architectures that work with this project
      systems = [
        "x86_64-linux" "aarch64-linux" "aarch64-darwin" "x86_64-darwin"
      ];
      perSystem = { config, self', inputs', pkgs, system, ... }: {

        # devShells.default describes the default shell with C++, cmake, boost,
        # and catch2
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            # C++ Compiler is already part of stdenv
            gnumake
            pkg-config
            boost
            cgal_5
            fftw
            gsl
            hdf5-cpp
            yaml-cpp
            fmt
            tbb
            gmp
            mpfr
          ];
        };
      };
    };
}
