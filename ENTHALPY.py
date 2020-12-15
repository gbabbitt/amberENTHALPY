#!/usr/bin/env python
import os
from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty

print("Welcome to AMBER BINDING ENTHALPY CALCULATOR- a simple graphical user interface for binding enthalpy estimation based upon GPU accelerated AMBER molecular dynamic simulations")
cmd = 'gedit READMEenthalpy.md'
os.system(cmd)
print("finding paths for paths.ctl") 
cmd = 'perl PATHSamber.pl'
os.system(cmd)


class ENTHALPYApp(App):
#    kv_directory = 'kivy_templates'
    def build(self):
        return MyLayout()
   
class MyLayout(Widget):
    
      
    # define buttons and actions
    def btn0(self):
            print("dual GPU mode selected")
            self.gpu = 2
            return self.gpu
    def btn00(self):
            print("single GPU mode selected")
            self.gpu = 1
            return self.gpu    
    def btn1(self):
        if self.gpu == 1:
            print("running amberENTHALPY on single GPU") 
            cmd = 'perl GUI_START_ENTHALPYlp1.pl'
            os.system(cmd)
        if self.gpu == 2:
            print("running amberENTHALPY on dual GPU") 
            cmd = 'perl GUI_START_ENTHALPYlp1_dualGPU.pl'
            os.system(cmd)
    
    


if __name__ == '__main__':
    ENTHALPYApp().run()
