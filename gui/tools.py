# Author: Lucas Zenichi Terada
# Institution: University of Campinas
# Description: This files contains functions to style

# Functions for fonts
def huge():
    return 'Verdana', 20, 'bold'


def large():
    return 'Verdana', 18


def normal():
    return 'Verdana', 12


# Default function for white background labels
def title():
    return {'font': huge(), 'background': 'white'}


def label():
    return {'font': large(), 'background': 'white'}


# Functions for buttons
def bgreen(text):
    return {'text': text, 'background': '#00e68a', 'font': large(), 'width': 15,
            'foreground': 'white'}


def bred(text):
    return {'text': text, 'background': '#ff0000', 'font': large(), 'width': 15,
            'foreground': 'white'}


def blred(text):
    return {'text': text, 'background': '#ff0000', 'font': normal(), 'width': 15,
            'foreground': 'white'}


def blblue(text):
    return {'text': text, 'background': '#0000ff', 'font': normal(), 'width': 15,
            'foreground': 'white'}


def blgreen(text):
    return {'text': text, 'background': '#00e68a', 'font': normal(), 'width': 15,
            'foreground': 'white'}
