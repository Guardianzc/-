import speech_recognition as sr
r =  sr.Recognizer()
with sr.WavFile(path) as source:
    audio = r.record(source)
IBM_USERNAME = "61c76251-869d-4911-9700-d07acdf178fe"
IBM_PASSWORD = "6AjpigA5qX5j"

#{
#  "url": "https://stream.watsonplatform.net/speech-to-text/api",
#  "username": "61c76251-869d-4911-9700-d07acdf178fe",
#  "password": "6AjpigA5qX5j"
#}

text = r.recognize_ibm(audio, username = IBM_USERNAME, password = IBM_PASSWORD, language = 'fr-FR') 