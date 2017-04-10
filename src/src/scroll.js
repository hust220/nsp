  function elementPosition(obj) {
    var curleft = 0
    var curtop = 0
    if (obj.offsetParent) {
      curleft = obj.offsetLeft
      curtop = obj.offsetTop
      obj = obj.offsetParent
      while (obj) {
        curleft += obj.offsetLeft
        curtop += obj.offsetTop
        obj = obj.offsetParent
      }
    }
    return { x: curleft, y: curtop }
  }

  export function ScrollToControl(id) {
    var h = document.documentElement.scrollTop
    console.log('id: ' + id)
    if (id === '') {
      window.scrollBy(0, 0 - h - 60)
    } else {
      var elem = document.getElementById(id)
      if (elem) {
        var scrollPos = elementPosition(elem).y
        scrollPos = scrollPos - h - 60
        console.log('document: ' + h + ', ele: ' + elementPosition(elem).y)
        window.scrollBy(0, scrollPos)
      }
    }
  }

  // var repeatCount = 0
  // var cTimeout
  // var timeoutIntervals = []
  // var timeoutIntervalSpeed

  // function ScrollSmoothly(scrollPos, repeatTimes) {
  //   if (repeatCount < repeatTimes) {
  //     window.scrollBy(0, 50)
  //   } else {
  //     repeatCount = 0
  //     clearTimeout(cTimeout)
  //     return
  //   }
  //   repeatCount++
  //   cTimeout = setTimeout(ScrollSmoothly(scrollPos, repeatTimes), 10)
  // }

