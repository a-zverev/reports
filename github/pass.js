var pass = prompt("Enter the Password:", "");
if (pass == null)  window.location = "fail.html";
else if (pass.toLowerCase() == "root")  
      window.location = "ok.html";
else  window.location = "fail.html";
