<template>
  <el-row class="doc">
    <el-col :span="6" class="sidebar">
      <div class="doc-nav" v-html="doc_nav"></div>
    </el-col>
    <el-col :span="18" :offset="6" class="main">
      <div class="doc-content" v-html="content"></div>
    </el-col>
  </el-row>
</template>

<script>
  import { bus } from '../bus.js'
  import marked from 'marked'

  function elementPosition(obj) {
    var curleft = 0
    var curtop = 0
    console.log(obj)
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

  function ScrollToControl(id) {
    console.log(id)
    var elem = document.getElementById(id)
    if (elem) {
      console.log(elem)
      var scrollPos = elementPosition(elem).y
      console.log(scrollPos)
      console.log(document.body.scrollTop)
      scrollPos = scrollPos - document.body.scrollTop - 60
      console.log(scrollPos)
      window.scrollBy(0, scrollPos)
      // var remainder = scrollPos % 50
      // var repeatTimes = (scrollPos - remainder) / 50
      // ScrollSmoothly(scrollPos, repeatTimes)
      // window.scrollBy(0, remainder)
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

  export default {
    data() {
      return {
        doc_nav: '',
        content: 'NSP'
      }
    },

    methods: {
      scrollTo(id) {
        console.log(id)
        if (id) {
          ScrollToControl(id)
        }
      },

      fetch_doc_nav() {
        var url = '/static/docs/nav.md'
        this.$http.get(url).then(response => {
          this.doc_nav = marked(response.body)
        }, response => {
          //
        })
      },

      fetch_content(theme, item) {
        var url = '/static/docs/' + theme + '.md'
        this.$http.get(url).then(response => {
          this.content = marked(response.body)
          this.scrollTo(item)
        }, response => {
          //
        })
      }
    },

    created() {
      this.fetch_doc_nav()
      var theme = (this.$route.params.theme ? this.$route.params.theme : 'install')
      var item = (this.$route.params.item ? this.$route.params.item : '')
      this.fetch_content(theme, item)
    },

    watch: {
      '$route'(to, from) {
        var theme = (to.params.theme ? to.params.theme : 'install')
        var item = (this.$route.params.item ? this.$route.params.item : '')
        this.fetch_content(theme, item)
      }
    },

  }
</script>

<style>
  .sidebar {
    position: fixed;
    font-size: 0.85em;
    overflow-y:scroll;
    height: 100%;
    padding-left: 50px;
    padding-bottom: 100px;
  }

  .sidebar ul {
    margin: 10px 15px;
    padding: 0px;
  }

  .sidebar li {
    list-style-type: none;
    margin: 10px 0px;
  }

  .doc-nav >ul {
    margin: 10px 0px;
  }

  .doc-content {
    margin: 35px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";
  }

  .doc-content pre {
    background-color: #f7f7f7;
    padding: 16px;
  }

  .doc-content code {
    font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  }

  .doc-content table {
    border-collapse: collapse;
  }

  .doc-content tr:nth-child(2n) {
    background-color: #f8f8f8;
  }

  .doc-content td {
    padding: 6px 13px;
    border: 1px solid #ddd;
  }

</style>