"""Module for displaying a menu to choose a node in a given hdf5 file"""
#import subprocess
import sys
import os
import signal
import urwid
import h5py


## nicer exit from ctrl-c
def signalHandler(signal, frame):
        print(' Keyboard interupt, exiting.')
        sys.exit(0)
signal.signal(signal.SIGINT, signalHandler)



PALETTE = [
    (None, 'light gray', 'black'),
    ('titlebar', 'white, bold', 'black'),
    ('footer', 'white', 'black'),
    ('heading', 'black', 'light gray'),
    ('line', 'black', 'light gray'),
    ('options', 'dark gray', 'black'),
    ('focus heading', 'white, bold', 'dark red'),
    ('focus line', 'black', 'dark red'),
    ('focus options', 'black', 'light gray'),
    ('selected', 'white', 'dark gray'),
    ('enter', 'dark green, bold', 'black'), 
    ('quit', 'dark red, bold', 'black')]

FOCUS_MAP = {
    'heading': 'focus heading',
    'options': 'focus options',
    'line': 'focus line'}


# {{{ the urwid classes


class MenuButton(urwid.Button):
    def __init__(self, caption, callback):
        super(MenuButton, self).__init__("")
        urwid.connect_signal(self, 'click', callback)
        self._w = urwid.AttrMap(urwid.SelectableIcon(
            [u'  \N{BULLET} ', caption], 2), None, 'selected')


class SubMenu(urwid.WidgetWrap):
    def __init__(self, path, caption, choices, top):
        super(SubMenu, self).__init__(MenuButton(
            [caption, u"\N{HORIZONTAL ELLIPSIS}"], self.open_menu))
        self.path = path
        self.top = top
        line = urwid.Divider(u'\N{LOWER ONE QUARTER BLOCK}')
        listbox = urwid.ListBox(urwid.SimpleFocusListWalker([
            urwid.AttrMap(urwid.Text([u"\n  ", caption.split("/")[-1]]), 'heading'), # put end of filename in heading
            urwid.AttrMap(line, 'line'),
            urwid.Divider()] + choices + [urwid.Divider()]))
        self.menu = urwid.AttrMap(listbox, 'options')

    def open_menu(self, _):
        self.top.open_box(self.menu)


class Choice(urwid.WidgetWrap):
    def __init__(self, path, caption, top):
        super(Choice, self).__init__(
            MenuButton(caption, self.item_chosen))
        self.path = path
        self.caption = caption
        self.top = top

    def item_chosen(self, _):
        response_txt = urwid.Text(['  Dataset:\n\n  ', self.caption, "\n"])
        response = urwid.AttrMap(response_txt, "titlebar")

        line = urwid.Divider(u'\N{LOWER ONE QUARTER BLOCK}')

        done_button = MenuButton(u'Confirm choice', self.set_chosen)
        v_padding = urwid.Padding(done_button, left = 1, right = 1)
        done = urwid.LineBox(v_padding)

        response_box = urwid.ListBox(urwid.SimpleFocusListWalker(
                                    [urwid.AttrMap(urwid.Text([u"\n  ", self.caption]), 'heading'), \
                                     urwid.AttrMap(line, 'line'), urwid.Divider()] + [done] + [urwid.Divider()
                                    ]))
        self.top.open_box(urwid.AttrMap(response_box, 'options'))

    def set_chosen(self, _):
        self.top.chosen_element = self.path
        raise urwid.ExitMainLoop()


class HorizontalBoxes(urwid.Columns):
    def __init__(self):
        super(HorizontalBoxes, self).__init__([], dividechars=1)
        self.chosen_element = None

    def open_box(self, box):
        if self.contents:
            del self.contents[self.focus_position + 1:]
        self.contents.append((urwid.AttrMap(box, 'options', FOCUS_MAP),
                              self.options('given', 24)))
        self.focus_position = len(self.contents) - 1
# }}}


# {{{ scanHDF5 and construct the tree as we go
def scan_hdf5(filename, top, tab_step=2):
    """Display the contents of a HDF5 file.

       Modified from: https://stackoverflow.com/questions/43371438/how-to-inspect-h5-file-in-python
    """

    def scan_node(node, tabs=0):
        path = node.name
        node_name = os.path.basename(path)
        if not node_name:
            node_name = node.filename

        if isinstance(node, h5py.Group):
            menu_list = []
            for entry in node.values():
                scan_result = scan_node(entry, tabs + tab_step)
                if scan_result is not None:
                    menu_list.append(scan_result)
            return SubMenu(path, node_name, menu_list, top)

        if isinstance(node, h5py.Dataset):
            return Choice(path, node_name, top)

        return None

    with h5py.File(filename, 'r') as hdf5_file:
        menu = scan_node(hdf5_file)
        return menu
# }}}


def handle_input(key):
    if key == 'Q' or key == 'q':
        raise urwid.ExitMainLoop()


def show_menu(filename):
    """Main entry point for the program"""

    # create a layout to pass to urwid.MainLoop()
    # -------------------------------------------

    # header
    #header_txt = urwid.Text([('', " Browsing file "), filename, ('', "\n Select a "), "vertex" ,("", " attribute to plot.")])
    header_txt = urwid.Text([" ", filename]) # , ('', "\n Select a "), "vertex" ,("", " attribute to plot.")])
    header = urwid.AttrMap(header_txt, "titlebar")

    # footer
    footer_txt = urwid.Text([u" Press ", ('enter', u'ENTER'), u' to select option.', u" Press ", ('quit', u'Q'), u' to quit.'])
    footer = urwid.AttrMap(footer_txt, "footer")

    # main body i.e. the menus
    top = HorizontalBoxes()
    menu_top = scan_hdf5(filename, top)
    top.open_box(menu_top.menu)

    menu = urwid.Filler(top, "middle", 10)
    v_padding = urwid.Padding(menu, left = 1, right = 1)
    box = urwid.LineBox(v_padding)

    #layout = urwid.Frame(header = header, body = menu)
    layout = urwid.Frame(header = header, body = box, footer = footer)

    # call the loop
    urwid.MainLoop(layout, PALETTE, unhandled_input = handle_input).run()

    #subprocess.call("clear", shell=True)
    print("Chosen:", top.chosen_element)

    return top.chosen_element

