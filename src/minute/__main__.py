#!/usr/bin/env python3
"""
Minute
"""
import ast
import logging
import pkgutil
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import importlib.resources

from . import cli, __version__
from .cli import CommandLineError


logger = logging.getLogger(__name__)


def main(arguments=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    subcommand_name = get_subcommand_name(arguments)
    module = importlib.import_module("." + subcommand_name, cli.__name__)

    parser = ArgumentParser(description=__doc__, prog="minute")
    parser.add_argument("--debug", action="store_true", help="Print some debugging information")
#    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser(
        subcommand_name, help=module.__doc__.split("\n")[1], description=module.__doc__, formatter_class=RawDescriptionHelpFormatter
    )
    module.add_arguments(subparser)
    args, remainder = parser.parse_known_args(arguments)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug
    try:
        module.main(args, arguments=remainder)
    except CommandLineError as e:
        logger.error(e)
        sys.exit(1)


def get_subcommand_name(arguments) -> str:
    """
    Parse arguments to find out which subcommand was requested.

    This sets up a minimal ArgumentParser with the correct help strings.

    Because help is obtained from a module’s docstring, but importing each module
    makes startup slow, the modules are only parsed with the ast module, and
    not fully imported at this stage.

    Return:
        subcommand name
    """
    parser = ArgumentParser(description=__doc__, prog="minute")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()

    for module_name, docstring in cli_modules():
        subparser = subparsers.add_parser(
            module_name, help=docstring.split("\n")[1], description=docstring, add_help=False
        )
        subparser.set_defaults(module_name=module_name)
    args, _ = parser.parse_known_args(arguments)
    module_name = getattr(args, "module_name", None)
    if module_name is None:
        parser.error("Please provide the name of a subcommand to run")
    return module_name


def cli_modules():
    """
    Yield (module_name, docstring) tuples for all modules in the "minutee.cli" package.
    """
    modules = pkgutil.iter_modules(cli.__path__)
    for module in modules:
        spec = importlib.util.find_spec(cli.__name__ + "." + module.name)
        with open(spec.origin) as f:
            mod_ast = ast.parse(f.read())
        docstring = ast.get_docstring(mod_ast, clean=False)
        yield module.name, docstring


if __name__ == "__main__":
    main()
