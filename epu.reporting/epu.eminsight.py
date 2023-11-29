#!/dls_sw/apps/EM/ebic.scripts/0.0/ebic.scripts_python/bin/python

import sys
import os
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, QLineEdit, QLabel, QTabWidget, QCheckBox
from PyQt5.QtCore import Qt, QProcess
import subprocess
import xmltodict
import json

class FileCheckerWidget(QWidget):
    def __init__(self, filename, is_directory=False):
        super().__init__()
        self.filename = filename
        self.is_directory = is_directory

        self.layout = QVBoxLayout()

        self.label = QLabel(self.filename)
        self.layout.addWidget(self.label, alignment=Qt.AlignLeft)

        self.status_label = QLabel()
        self.layout.addWidget(self.status_label, alignment=Qt.AlignRight)

        self.path_label = QLabel()
        self.layout.addWidget(self.path_label)

        self.setLayout(self.layout)

    def check_directory(self, selected_script):
        executable_path = os.path.dirname(selected_script)
        full_path = os.path.join(executable_path, self.filename)
        return os.path.isdir(full_path), full_path

    def check_file(self, selected_script):
        executable_path = os.path.dirname(selected_script)
        full_path = os.path.join(executable_path, self.filename)
        return os.path.exists(full_path), full_path

    def readXmlSession(self, xml_path):
        # It is necessary to have a function for getting xml session name elsewhere in script
        with open(xml_path, "r") as xml:
            for_parsing = xml.read()
            data = xmltodict.parse(for_parsing)
        data = data["EpuSessionXml"]

        sessionName = data["Name"]["#text"]
        
        return sessionName
    
    def readAtlasSession(self, xml_path):
        # This will fetch the atlas xml data
        with open(xml_path, "r") as xml:
            for_parsing = xml.read()
            data = xmltodict.parse(for_parsing)
        data = data["ScreeningSessionXml"]

        sessionName = data["Name"]['#text']

        return sessionName

    def update_file_status(self, selected_script):
        if self.is_directory:
            exists, path = self.check_directory(selected_script)
        else:
            exists, path = self.check_file(selected_script)

        if exists:
            self.status_label.setText("Found")
            self.status_label.setStyleSheet("QLabel { color: green; }")
            self.path_label.setText("Path: " + path)
        else:
            self.status_label.setText("Not Found")
            self.status_label.setStyleSheet("QLabel { color: red; }")
            self.path_label.setText("")

    def update_file_status_session(self, selected_script):
        if self.is_directory:
            exists, path = self.check_directory(selected_script)
        else:
            exists, path = self.check_file(selected_script)

        if exists:
            sessionName = self.readXmlSession(selected_script)
            #sessionName = ":"
            self.status_label.setText("Found: "+str(sessionName))
            self.status_label.setStyleSheet("QLabel { color: green; }")
            self.path_label.setText("Path: " + path)
        else:
            self.status_label.setText("Not Found")
            self.status_label.setStyleSheet("QLabel { color: red; }")
            self.path_label.setText("")


    def update_file_status_atlas(self, selected_script):
        if self.is_directory:
            exists, path = self.check_directory(selected_script)
        else:
            exists, path = self.check_file(selected_script)

        if exists:
            sessionName = self.readAtlasSession(selected_script)
            #sessionName = ":"
            self.status_label.setText("Found: "+str(sessionName))
            self.status_label.setStyleSheet("QLabel { color: green; }")
            self.path_label.setText("Path: " + path)
        else:
            self.status_label.setText("Not Found")
            self.status_label.setStyleSheet("QLabel { color: red; }")
            self.path_label.setText("")
   

class SessionTab(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EMinsight")
        self.resize(570, 310)  # Set the default window size here

        self.dark_mode = True  # Set dark mode as the default appearance
        self.light_mode_stylesheet = ""
        self.dark_mode_stylesheet = "background-color: #333; color: #fff;"

        self.layout = QVBoxLayout()

        self.file_button = QPushButton("Select EPU session: EpuSession.dm")
        self.file_button.clicked.connect(self.open_file_dialog)
        self.layout.addWidget(self.file_button)

        self.script_path = QLineEdit()
        self.layout.addWidget(self.script_path)

        self.file_layout = QHBoxLayout()
        self.file_layout.addWidget(self.file_button)
        self.file_layout.addWidget(self.script_path)

        self.file_layout.setStretch(0, 1)
        self.file_layout.setStretch(1, 0)

        self.layout.addLayout(self.file_layout)

        self.file_checkers_layout = QVBoxLayout()

        self.python_script_checker = FileCheckerWidget("EpuSession.dm", is_directory=False)
        self.file_checkers_layout.addWidget(self.python_script_checker)

        self.data_file_checker = FileCheckerWidget("Images-Disc1", is_directory=True)
        self.file_checkers_layout.addWidget(self.data_file_checker)

        # Disable the third file checker by not adding it to the layout
        # self.config_file_checker = FileCheckerWidget("ScreeningSession.dm", is_directory=False)
        # self.file_checkers_layout.addWidget(self.config_file_checker)

        self.layout.addLayout(self.file_checkers_layout)

        # Adding the second file browser button underneath the first one
        self.second_file_button = QPushButton("Select Atlas session: ScreeningSession.dm")
        self.second_file_button.clicked.connect(self.open_second_file_dialog)
        self.layout.addWidget(self.second_file_button)

        self.second_file_path = QLineEdit()
        self.layout.addWidget(self.second_file_path)

        # Add a fourth FileCheckerWidget for the file selected using the "Select Another File" button
        self.fourth_file_checker = FileCheckerWidget("ScreeningSession.dm", is_directory=False)
        self.layout.addWidget(self.fourth_file_checker)

        self.fourth_file_path = QLabel()
        self.layout.addWidget(self.fourth_file_path)

        # Remove the third FileCheckerWidget and its label from the layout
        # self.third_file_checkers_layout = QVBoxLayout()
        # self.third_file_checker = FileCheckerWidget("ScreeningSession.dm", is_directory=False)
        # self.third_file_checkers_layout.addWidget(self.third_file_checker)

        # self.third_file_path = QLabel()
        # self.third_file_checkers_layout.addWidget(self.third_file_path)

        # self.layout.addLayout(self.third_file_checkers_layout)

        # Add a label and text box for command line options
        self.command_line_label = QLabel("Command line options:")
        self.layout.addWidget(self.command_line_label)

        # Command line options
        self.command_line_textbox = QLineEdit()
        self.layout.addWidget(self.command_line_textbox)

        # Square analysis checkbox
        self.gridsquareanalysis_checkbox = QCheckBox("GridSquare analysis")
        self.layout.addWidget(self.gridsquareanalysis_checkbox)
        
        # Development mode analysis checkbox
        self.devmode_checkbox = QCheckBox("dev mode")
        self.layout.addWidget(self.devmode_checkbox)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.gridsquareanalysis_checkbox)
        button_layout.addWidget(self.devmode_checkbox)

        # Adding the Deposition button next to the Execute button
        self.deposition_button = QPushButton("Prepare deposition file")
        self.deposition_button.clicked.connect(self.deposition_script)

        # Adding the execute button
        self.execute_button = QPushButton("Execute EMinsight")
        self.execute_button.clicked.connect(self.execute_script)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.deposition_button)
        button_layout.addWidget(self.execute_button)

        self.layout.addLayout(button_layout)

        # Save settings buttons
        self.save_button = QPushButton("Save Settings")
        self.load_button = QPushButton("Load Settings")

        self.save_button.clicked.connect(self.save_settings)
        self.load_button.clicked.connect(self.load_settings)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.load_button)

        self.layout.addLayout(button_layout)

        self.setLayout(self.layout)

        # Set the initial stylesheet to dark mode
        self.setStyleSheet(self.dark_mode_stylesheet)

    def save_settings(self):
        settings = {
            "data_directory": self.script_path.text(),
            "output_directory": self.second_file_path.text(),
            "gridsquare_analysis_checkbox": self.gridsquareanalysis_checkbox.isChecked(),
            "devmode_checkbox": self.devmode_checkbox.isChecked()
        }

        file_dialog = QFileDialog()
        options = file_dialog.Options()
        options |= file_dialog.DontUseNativeDialog
        default_filename = "eminsight_settings_session.json"
        file_name, _ = file_dialog.getSaveFileName(self, "Save Settings to File", default_filename, "Settings File (*.json);;All Files (*)", options=options)

        if file_name:
            with open(file_name, 'w') as file:
                json.dump(settings, file)

    def load_settings(self):
        file_dialog = QFileDialog()
        options = file_dialog.Options()
        options |= file_dialog.DontUseNativeDialog
        file_name, _ = file_dialog.getOpenFileName(self, "Load Settings from File", os.getcwd(), "Settings File (*.json);;All Files (*)", options=options)

        if file_name:
            try:
                with open(file_name, 'r') as file:
                    settings = json.load(file)

                    self.script_path.setText(settings.get("data_directory", ""))
                    self.second_file_path.setText(settings.get("output_directory", ""))
                    self.gridsquareanalysis_checkbox.setChecked(settings.get("gridsquare_analysis_checkbox", False))
                    self.devmode_checkbox.setChecked(settings.get("devmode_checkbox", False))
                
                self.selected_script = self.script_path.text()
                self.update_file_status_session()
                self.second_selected_script = self.second_file_path.text()
                self.update_file_status_atlas()
            except Exception as e:
                print("Error loading settings:", e)

    def toggle_dark_mode(self):
        if self.dark_mode:
            self.setStyleSheet(self.light_mode_stylesheet)
        else:
            self.setStyleSheet(self.dark_mode_stylesheet)
        self.dark_mode = not self.dark_mode

    def open_file_dialog(self):
        initial_dir = os.getcwd()
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        filename, _ = file_dialog.getOpenFileName(self, "Select EPU session: EpuSession.dm", initial_dir)
        if filename:
            self.selected_script = filename
            self.script_path.setText(filename)
            self.update_file_status_session()

    def open_second_file_dialog(self):
        initial_dir = os.getcwd()
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        filename, _ = file_dialog.getOpenFileName(self, "Select Atlas session: ScreeningSession.dm", initial_dir)
        if filename:
            self.second_selected_script = filename
            self.second_file_path.setText(filename)
            self.update_file_status_atlas()

    def update_file_status_session(self):
        self.python_script_checker.update_file_status_session(self.selected_script)
        self.data_file_checker.update_file_status_session(self.selected_script)
        # Skip the update for the third file checker since it is disabled
        # self.config_file_checker.update_file_status(self.selected_script)

    def update_file_status_atlas(self):
        self.fourth_file_checker.update_file_status_atlas(self.second_selected_script)

    def update_file_statuses(self):
        self.python_script_checker.update_file_status(self.selected_script)
        self.data_file_checker.update_file_status(self.selected_script)
        # Skip the update for the third file checker since it is disabled
        # self.config_file_checker.update_file_status(self.selected_script)

    def execute_script(self):
        print('Running full EMinsight analysis')
        try:
            # Check if dev script needs to be run
            devmode = self.devmode_checkbox.isChecked()

            script_dir = os.path.dirname(os.path.abspath(__file__))
            if devmode:
                script_to_execute = os.path.join(script_dir, "epu.xml_summary_dev.py")
            else:
                script_to_execute = os.path.join(script_dir, "epu.xml_summary.py")

            # Pull values from text box GUI in case they have been modified
            self.selected_script = self.script_path.text()
            self.second_selected_script = self.second_file_path.text()
            command_line_options = self.command_line_textbox.text()

            # Check if reanalyse checkbox is checked
            gridsquareanalysis = self.gridsquareanalysis_checkbox.isChecked()

            # Modify the command_line_args based on the gridsquareanalysis checkbox
            if gridsquareanalysis:
                command_line_args = command_line_options + ' --screening Y'
            else:
                command_line_args = command_line_options

            if os.path.exists(script_to_execute):
                command_line_args = '--i ' + str(self.selected_script) + ' --o processing --atlas ' + str(
                    self.second_selected_script) + ' ' + str(command_line_args)

                # Use subprocess.call() to execute the script with command line arguments
                print('+++ '+str(script_to_execute) + ' ' + str(command_line_args))
                subprocess.call(str(script_to_execute) + ' ' + str(command_line_args), shell=True)
            else:
                print("Script not found:", script_to_execute)
        except AttributeError:
            print("Please select EpuSession.dm and Atlas ScreeningSession.dm first.")

    def deposition_script(self):
        print('Preparing deposition file')
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            script_to_execute = os.path.join(script_dir, "epu.xml_summary.py")

            # Pull values from text box GUI in case they have been modified
            self.selected_script = self.script_path.text()
            self.second_selected_script = self.second_file_path.text()
            command_line_options = self.command_line_textbox.text()

            # Force deposition flag
            command_line_args = '--doppio Y'

            if os.path.exists(script_to_execute):
                # Replace the command line arguments as needed for the deposition script
                command_line_args = '--i ' + str(self.selected_script) + ' --o processing --atlas ' + str(
                    self.second_selected_script) + ' ' + str(command_line_args)

                # Use subprocess.call() to execute the script with command line arguments
                print('+++ '+str(script_to_execute) + ' ' + str(command_line_args))
                subprocess.call(str(script_to_execute) + ' ' + str(command_line_args), shell=True)
            else:
                print("Script not found:", script_to_execute)
        except AttributeError:
            print("Please select EpuSession.dm and Atlas ScreeningSession.dm first.")

class GlobalTab(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()

        # Add a label and text box for the data directory at the top
        self.select_data_button = QPushButton("Select data directory")
        self.select_data_button.clicked.connect(self.open_data_directory_dialog)
        self.layout.addWidget(self.select_data_button)

        self.data_directory_textbox = QLineEdit()
        self.layout.addWidget(self.data_directory_textbox)

        # Add a label and text box for the output directory at the top
        self.select_output_button = QPushButton("Select output directory")
        self.select_output_button.clicked.connect(self.open_output_directory_dialog)
        self.layout.addWidget(self.select_output_button)

        self.output_directory_textbox = QLineEdit()
        self.layout.addWidget(self.output_directory_textbox)

        # Add a label and text box for command line options
        self.command_line_label = QLabel("Command line options:")
        self.layout.addWidget(self.command_line_label)

        self.command_line_textbox = QLineEdit()
        self.layout.addWidget(self.command_line_textbox)

        # Save settings buttons
        self.save_button = QPushButton("Save Settings")
        self.load_button = QPushButton("Load Settings")

        self.save_button.clicked.connect(self.save_settings)
        self.load_button.clicked.connect(self.load_settings)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.load_button)

        # Reanalyse checkbox
        self.reanalyse_checkbox = QCheckBox("reanalyse existing")
        self.layout.addWidget(self.reanalyse_checkbox)

        # Add a button for executing the external script
        self.execute_button = QPushButton("Execute EMinsight global analysis")
        self.execute_button.clicked.connect(self.execute_external_script)
        self.layout.addWidget(self.execute_button)

        self.layout.addLayout(button_layout)
        self.setLayout(self.layout)

    def open_data_directory_dialog(self):
        initial_dir = os.getcwd()
        directory_path = QFileDialog.getExistingDirectory(self, "Select Output Directory", initial_dir)
        if directory_path:
            self.data_directory_textbox.setText(directory_path)
            self.selected_dir_in = directory_path

    def open_output_directory_dialog(self):
        initial_dir = os.getcwd()
        directory_path = QFileDialog.getExistingDirectory(self, "Select Output Directory", initial_dir)
        if directory_path:
            self.output_directory_textbox.setText(directory_path)
            self.selected_dir_out = directory_path

    def execute_external_script(self):
        print('Running global analysis')
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            script_to_execute = os.path.join(script_dir, "epu.reporting.py")

            # Pull values from text box GUI in case they have been modified
            self.selected_dir_in = self.data_directory_textbox.text()
            self.selected_dir_out = self.output_directory_textbox.text()
            command_line_options = self.command_line_textbox.text()
            
            # Check if reanalyse checkbox is checked
            reanalyse = self.reanalyse_checkbox.isChecked()

            # Modify the command_line_args based on the reanalyse checkbox
            if reanalyse:
                command_line_args = command_line_options + ' --reanalyse Y'
            else:
                command_line_args = command_line_options

            if os.path.exists(script_to_execute):
                # Retrive any user input command line options
                command_line_args = '--input ' + str(self.selected_dir_in) + ' --output ' + str(self.selected_dir_out) + ' ' + str(command_line_args)

                # Use subprocess.call() to execute the script with command line arguments
                print('+++ '+str(script_to_execute) + ' ' + str(command_line_args))
                subprocess.call(str(script_to_execute) + ' ' + str(command_line_args), shell=True)
            else:
                print("Script not found:", script_to_execute)
        except AttributeError:
            print("Please select an input directory and Atlas ScreeningSession.dm first.")

    def save_settings(self):
        settings = {
            "data_directory": self.data_directory_textbox.text(),
            "output_directory": self.output_directory_textbox.text(),
            "command_line_options": self.command_line_textbox.text(),
        }

        file_dialog = QFileDialog()
        options = file_dialog.Options()
        options |= file_dialog.DontUseNativeDialog
        default_filename = "eminsight_settings_global.json"
        file_name, _ = file_dialog.getSaveFileName(self, "Save Settings to File", default_filename, "Settings File (*.json);;All Files (*)", options=options)

        if file_name:
            with open(file_name, 'w') as file:
                json.dump(settings, file)

    def load_settings(self):
        file_dialog = QFileDialog()
        options = file_dialog.Options()
        options |= file_dialog.DontUseNativeDialog
        file_name, _ = file_dialog.getOpenFileName(self, "Load Settings from File", os.getcwd(), "Settings File (*.json);;All Files (*)", options=options)

        if file_name:
            try:
                with open(file_name, 'r') as file:
                    settings = json.load(file)

                    self.data_directory_textbox.setText(settings.get("data_directory", ""))
                    self.output_directory_textbox.setText(settings.get("output_directory", ""))
                    self.command_line_textbox.setText(settings.get("command_line_options", ""))

            except Exception as e:
                print("Error loading settings:", e)

class MainWidget(QWidget):
    def __init__(self):
        super().__init__()

        self.dark_mode = True  # Set dark mode as the default appearance
        self.light_mode_stylesheet = ""
        self.dark_mode_stylesheet = "background-color: #333; color: #fff;"

        # Adding the dark mode toggle button to the top of the layout
        self.toggle_mode_button = QPushButton("Toggle Dark Mode")
        self.toggle_mode_button.clicked.connect(self.toggle_dark_mode)

        self.layout = QVBoxLayout()
        self.layout.addWidget(self.toggle_mode_button, alignment=Qt.AlignRight)  # Move the button to the top

        # Create the QTabWidget
        self.tabs = QTabWidget()

        # Create the 'Session' and 'Global' tabs
        self.session_tab = SessionTab()
        self.global_tab = GlobalTab()

        # Add the 'Session' and 'Global' tabs to the QTabWidget
        self.tabs.addTab(self.session_tab, "Session")
        self.tabs.addTab(self.global_tab, "Global")

        # Set the dark mode stylesheet for the tab widget and each tab
        self.tabs.setStyleSheet("QTabWidget::pane { background-color: #333; }"
                                "QTabBar::tab { color: #fff; }"
                                "QTabBar::tab:selected { background-color: #444; }")
        self.session_tab.setStyleSheet(self.dark_mode_stylesheet)  # Set initial dark mode stylesheet

        # Add the tab widget to the main layout
        self.layout.addWidget(self.tabs)

        self.setLayout(self.layout)

        # Set the initial stylesheet to dark mode
        self.setStyleSheet(self.dark_mode_stylesheet)

    def toggle_dark_mode(self):
        self.dark_mode = not self.dark_mode
        if self.dark_mode:
            self.setStyleSheet(self.dark_mode_stylesheet)
            self.session_tab.setStyleSheet(self.dark_mode_stylesheet)
            self.global_tab.setStyleSheet(self.dark_mode_stylesheet)
        else:
            self.setStyleSheet(self.light_mode_stylesheet)
            self.session_tab.setStyleSheet(self.light_mode_stylesheet)
            self.global_tab.setStyleSheet(self.light_mode_stylesheet)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWidget()
    window.show()
    sys.exit(app.exec_())