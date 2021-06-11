#!/usr/bin/envpython3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Comand line handlers
"""
import click

from cli.ats import ats
from cli.isoform import isoform
from cli.preprocess import preprocess
from cli.postprocess import postprocess


VERSION = "0.0.1-alpha"
LABEL = "scATS"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    context_settings=CONTEXT_SETTINGS,
)
@click.version_option(VERSION, message="Current version %(version)s")
def cli():
    u"""
    Welcome

    \b
    This function is used to test the function of sashimi plotting

    \f
    Created by ygidtu@gmail.com at 2018.12.19
    :return:
    """

    pass


cli.add_command(preprocess)
cli.add_command(ats)
cli.add_command(isoform)
cli.add_command(postprocess)
