let
  pkgs = import <nixpkgs> {};
  stdenv = pkgs.overrideCC pkgs.stdenv pkgs.gcc7;
in rec {
  ecdlp-diem = stdenv.mkDerivation rec {
    name = "ecdlp-diem" ;
    version = "0.1" ;
    buildInputs =
      with pkgs; [
        cmake gnumake
      ] ;
  } ;
}

